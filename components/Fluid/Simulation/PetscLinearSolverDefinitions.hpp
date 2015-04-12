#include <Eigen/Core>

#include <petsc.h>

#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class VpeStencilGenerator {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridGeometryType = typename SolverTraitsType::GridGeometryType;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  VpeStencilGenerator() {}

  void
  initialize(ParallelDistributionType const* parallelDistribution,
             ParametersType const*           parameters,
             ScalarType const*               dt) {
    _parallelDistribution = parallelDistribution;
    _parameters           = parameters;
    _dt                   = dt;
  }

  template <typename TStencil,
            typename TRow,
            typename TColumns>
  inline void
  get(CellAccessorType const& accessor,
      TStencil&               stencil,
      TRow&                   row,
      TColumns&               columns) const {
    typedef Eigen::Matrix<ScalarType, 2* Dimensions, 1> Vector2Ds;

    auto corner = _parallelDistribution->corner;

    Vector2Ds meanWidths;

    for (int d = 0; d < Dimensions; ++d) {
      meanWidths(2 * d) = 0.5 *
                          (accessor.width(d) +
                           accessor.width(d, -1, d));
      meanWidths(2 * d + 1) = 0.5 *
                              (accessor.width(d) +
                               accessor.width(d, +1, d));

      auto leftIndex  = (accessor.relativeIndex(d, -1) + corner).eval();
      auto rightIndex = (accessor.relativeIndex(d, +1) + corner).eval();

      columns[2 * d].i     = leftIndex(0);
      columns[2 * d].j     = leftIndex(1);
      columns[2 * d + 1].i = rightIndex(0);
      columns[2 * d + 1].j = rightIndex(1);

      if (Dimensions == 3) {
        columns[2 * d].k     = leftIndex(2);
        columns[2 * d + 1].k = rightIndex(2);
      }
    }

    auto currentIndex = (accessor.index() + corner).eval();

    columns[2 * Dimensions].i = currentIndex(0);
    columns[2 * Dimensions].j = currentIndex(1);

    if (Dimensions == 3) {
      columns[2 * Dimensions].k = currentIndex(2);
    }

    row = columns[2 * Dimensions];

    stencil[2 * Dimensions] = 1.0;

    ScalarType const coeff =  -(*_dt) / (2.0 * _parameters->re());

    for (int d = 0; d < Dimensions; ++d) {
      auto meanLeftAndRightWidth = meanWidths(2 * d) + meanWidths(2 * d + 1);

      stencil[2 * d]
        = 2.0 * coeff / (meanWidths(2 * d) * meanLeftAndRightWidth);
      stencil[2 * d + 1]
        = 2.0 * coeff / (meanWidths(2 * d + 1) * meanLeftAndRightWidth);

      ScalarType diagonalCoeff
        = -2.0 / (meanWidths(2 * d) * meanWidths(2 * d + 1));

      stencil[2 * Dimensions] += coeff * diagonalCoeff;
    }
    // logInfo("{1} {2} {3} {4} {5} {6}",
    // stencil[0],
    // stencil[1],
    // stencil[2],
    // stencil[3],
    // stencil[4],
    // Dimensions);
  }

private:
  ParallelDistributionType const* _parallelDistribution;
  ParametersType const*           _parameters;
  ScalarType const*               _dt;
};

template <unsigned TDimension>
class VpeRhsGenerator {
public:
  template <typename TCellAccessor>
  inline typename TCellAccessor::ScalarType
  get(TCellAccessor const& accessor) const {
    typename TCellAccessor::ScalarType result = accessor.fgh(TDimension);

    return result;
  }
};

template <int TDimension>
class VpeResultAcquirer {
public:
  template <typename TCellAccessor>
  inline void
  set(TCellAccessor const&                      accessor,
      typename TCellAccessor::ScalarType const& value) const {
    accessor.fgh(TDimension) = value;
  }
};

template <typename TSolverTraits>
class PpeStencilGenerator {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  PpeStencilGenerator() {}

  void
  initialize(ParallelDistributionType const* parallelDistribution) {
    _parallelDistribution = parallelDistribution;
  }

  template <typename TStencil,
            typename TRow,
            typename TColumns>
  inline void
  get(CellAccessorType const& accessor,
      TStencil&               stencil,
      TRow&                   row,
      TColumns&               columns) const {
    typedef Eigen::Matrix<ScalarType, 2* Dimensions, 1> Vector2Ds;

    auto corner = _parallelDistribution->corner;

    Vector2Ds meanWidths;

    for (int d = 0; d < Dimensions; ++d) {
      meanWidths(2 * d) = 0.5 *
                          (accessor.width(d) +
                           accessor.width(d, -1, d));
      meanWidths(2 * d + 1) = 0.5 *
                              (accessor.width(d) +
                               accessor.width(d, +1, d));

      auto leftIndex  = (accessor.relativeIndex(d, -1) + corner).eval();
      auto rightIndex = (accessor.relativeIndex(d, +1) + corner).eval();

      columns[2 * d].i     = leftIndex(0);
      columns[2 * d].j     = leftIndex(1);
      columns[2 * d + 1].i = rightIndex(0);
      columns[2 * d + 1].j = rightIndex(1);

      if (Dimensions == 3) {
        columns[2 * d].k     = leftIndex(2);
        columns[2 * d + 1].k = rightIndex(2);
      }
    }

    auto currentIndex = (accessor.index() + corner).eval();
    columns[2 * Dimensions].i = currentIndex(0);
    columns[2 * Dimensions].j = currentIndex(1);

    if (Dimensions == 3) {
      columns[2 * Dimensions].k = currentIndex(2);
    }

    row = columns[2 * Dimensions];

    stencil[2 * Dimensions] = 0;

    for (int d = 0; d < Dimensions; ++d) {
      auto meanLeftAndRightWidth = meanWidths(2 * d) + meanWidths(2 * d + 1);
      stencil[2 * d]
        = 2.0 / (meanWidths(2 * d) * meanLeftAndRightWidth);
      stencil[2 * d + 1]
        = 2.0 / (meanWidths(2 * d + 1) * meanLeftAndRightWidth);
      stencil[2 * Dimensions]
        -= 2.0 / (meanWidths(2 * d) * meanWidths(2 * d + 1));
    }
  }

private:
  ParallelDistributionType const* _parallelDistribution;
};

template <typename TSolverTraits>
class PpeRhsGenerator {
public:
  using CellAccessorType = typename TSolverTraits::CellAccessorType;

  using ScalarType = typename TSolverTraits::ScalarType;

  void
  initialize(ScalarType const* dt) { _dt = dt; }

  inline ScalarType
  get(CellAccessorType const& accessor) const {
    ScalarType result = 0.0;

    for (int d = 0; d < CellAccessorType::Dimensions; ++d) {
      result += (accessor.fgh(d) - accessor.fgh(d, -1, d))
                / accessor.width(d);
    }

    result = result / (*_dt);

    // logInfo("Rhs | {1} | {2}", accessor.index().transpose(), result);

    return result;
  }

private:
  ScalarType const* _dt;
};

class PpeResultAcquirer1 {
public:
  template <typename TCellAccessor>
  inline void
  set(TCellAccessor const&                      accessor,
      typename TCellAccessor::ScalarType const& value) const {
    accessor.pressure() = value;
  }
};

template <typename TSolverTraits>
class PpeResultAcquirer2 {
public:
  using CellAccessorType = typename TSolverTraits::CellAccessorType;

  using ScalarType = typename TSolverTraits::ScalarType;

  void
  initialize(ScalarType const* dt) { _dt = dt; }

  inline void
  set(CellAccessorType const& accessor,
      ScalarType const&       value) const {
    accessor.pressure()      += value;
    accessor.projectionTerm() = value;
  }

private:
  ScalarType const* _dt;
};
}
}
