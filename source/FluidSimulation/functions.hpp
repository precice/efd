#pragma once

#include "Private/mpigenerics.hpp"

#include <Eigen/Core>

#include <Uni/Logging/macros>

#include <petsc.h>

#include <cmath>
#include <limits>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TCellAccessor, typename TScalar, int TD>
void
computeMaxVelocity(TCellAccessor const&                  accessor,
                   typename TCellAccessor::VectorDsType& maxVelocity) {
  maxVelocity
    = accessor.velocity().cwiseAbs().cwiseQuotient(
    accessor.width()).cwiseMax(maxVelocity);
}

template <typename TScalar, int TD>
struct TimeStepProcessing {
  typedef Eigen::Matrix<TScalar, TD, 1> VectorDs;

  static inline TScalar
  compute(TScalar const&  re,
          TScalar const&  tau,
          VectorDs const& minCellWidth,
          VectorDs const& maxVelocity) {
    TScalar localMin, globalMin;
    TScalar factor;

    factor = minCellWidth.cwiseProduct(minCellWidth).cwiseInverse().sum();

    localMin = std::min((TScalar)(re / (2.0 * factor)),
                        maxVelocity.cwiseInverse().minCoeff());

    globalMin = std::numeric_limits<TScalar>::max();
    Private::mpiAllReduce<TScalar>(&localMin,
                                   &globalMin,
                                   1,
                                   MPI_MIN,
                                   PETSC_COMM_WORLD);

    factor  = globalMin;
    factor *= tau;

    return factor;
  }
};

template <typename TScalar>
inline TScalar
dudx(TScalar const& leftU,
     TScalar const& rightU,
     TScalar const& leftX,
     TScalar const& rightX) {
  TScalar const width = leftX + rightX;

  return (rightU - leftU) / width;
}

template <typename TScalar>
inline TScalar
d2udx2(TScalar const& currentU,
       TScalar const& leftU,
       TScalar const& rightU,
       TScalar const& currentX,
       TScalar const& rightX) {
  TScalar const factor = currentX + rightX;

  return 2.0 * (rightU / (rightX * factor) -
                currentU / (currentX * rightX) +
                leftU / (currentX * factor));
}

template <typename TScalar>
inline TScalar
d2udy2(TScalar const& currentU,
       TScalar const& bottomU,
       TScalar const& topU,
       TScalar const& currentY,
       TScalar const& bottomY,
       TScalar const& topY) {
  TScalar const factorBottom = 0.5 * (currentY + bottomY);
  TScalar const factorTop    = 0.5 * (currentY + topY);
  TScalar const factor       = factorBottom + factorTop;

  return 2.0 * (topU / (factorTop * factor) -
                currentU / (factorBottom * factorTop) +
                bottomU / (factorBottom * factor));
}

template <typename TScalar>
inline TScalar
duvdy(TScalar const& currentU,
      TScalar const& currentV,
      TScalar const& bottomU,
      TScalar const& bottomV,
      TScalar const& rightBottomV,
      TScalar const& rightV,
      TScalar const& topU,
      TScalar const& currentY,
      TScalar const& currentX,
      TScalar const& bottomY,
      TScalar const& topY,
      TScalar const& rightX,
      TScalar const& gamma) {
  TScalar const yCurrent = 0.5 * currentY;
  TScalar const yBottom  = 0.5 * (currentY + bottomY);
  TScalar const yTop     = 0.5 * (currentY + topY);
  TScalar const xCurrent = 0.5 * currentX;
  TScalar const xRight   = 0.5 * (currentX + rightX);

  TScalar const secondOrder
    = (((xRight - xCurrent) / xRight * currentV + xCurrent / xRight * rightV) *
       ((yTop - yCurrent) / yTop * currentU + yCurrent / yTop * topU) -
       ((xRight - xCurrent) / xRight * bottomV + xCurrent / xRight *
        rightBottomV) *
       ((yBottom - yCurrent) / yBottom * currentU + yCurrent / yBottom *
        bottomU)
       ) / (2.0 * yCurrent);

  TScalar const kr = (xRight - xCurrent) / xRight * currentV +
                     xCurrent / xRight * rightV;
  TScalar const kl = (xRight - xCurrent) / xRight * bottomV +
                     xCurrent / xRight * rightBottomV;

  TScalar const firstOrder = 1.0 / (4.0 * yCurrent) * (
    kr * (currentU + topU) - kl * (bottomU + currentU) +
    std::fabs(kr) * (currentU - topU) -
    std::fabs(kl) * (bottomU - currentU));

  TScalar const result = (1.0 - gamma) * secondOrder + gamma * firstOrder;

  return result;
}

template <typename TScalar>
inline TScalar
du2dx(TScalar const& currentU,
      TScalar const& leftU,
      TScalar const& rightU,
      TScalar const& currentX,
      TScalar const& leftX,
      TScalar const& rightX,
      TScalar const& gamma) {
  TScalar const dxShort = 0.5 * currentX;
  TScalar const dxLong0 = 0.5 * (leftX + currentX);
  TScalar const dxLong1 = 0.5 * (currentX + rightX);

  TScalar const kr
    = (dxLong1 - dxShort) / dxLong1 * currentU + dxShort / dxLong1 * rightU;

  TScalar const kl
    = (dxLong0 - dxShort) / dxLong0 * currentU + dxShort / dxLong0 * leftU;

  TScalar const secondOrder = (kr * kr - kl * kl) / (2.0 * dxShort);

  TScalar const firstOrder = 1.0 / (4.0 * dxShort) *
                             (kr * (currentU + rightU) -
                              kl * (leftU + currentU) +
                              std::fabs(kr) * (currentU - rightU) -
                              std::fabs(kl) * (leftU - currentU));

  TScalar const result = (1.0 - gamma) * secondOrder + gamma * firstOrder;

  return result;
}

template <typename TScalar>
inline TScalar
computeDiffusion2D(TScalar const& currentU,
                   TScalar const& leftU,
                   TScalar const& rightU,
                   TScalar const& bottomU,
                   TScalar const& topU,
                   TScalar const& currentX,
                   TScalar const& currentY,
                   TScalar const& rightX,
                   TScalar const& bottomY,
                   TScalar const& topY) {
  return d2udx2(currentU,
                leftU,
                rightU,
                currentX,
                rightX)
         + d2udy2(currentU,
                  bottomU,
                  topU,
                  currentY,
                  bottomY,
                  topY);
}

template <typename TScalar>
inline TScalar
computeDiffusion3D(TScalar const& currentU,
                   TScalar const& leftU,
                   TScalar const& rightU,
                   TScalar const& bottomU,
                   TScalar const& topU,
                   TScalar const& backU,
                   TScalar const& frontU,
                   TScalar const& currentX,
                   TScalar const& currentY,
                   TScalar const& currentZ,
                   TScalar const& rightX,
                   TScalar const& bottomY,
                   TScalar const& topY,
                   TScalar const& backZ,
                   TScalar const& frontZ) {
  return d2udx2(currentU,
                leftU,
                rightU,
                currentX,
                rightX)
         + d2udy2(currentU,
                  bottomU,
                  topU,
                  currentY,
                  bottomY,
                  topY)
         + d2udy2(currentU, // d2udz2
                  backU,
                  frontU,
                  currentZ,
                  backZ,
                  frontZ);
}

template <typename TScalar>
inline TScalar
computeConvection2D(TScalar const& currentU,
                    TScalar const& currentV,
                    TScalar const& leftU,
                    TScalar const& rightU,
                    TScalar const& rightV,
                    TScalar const& bottomU,
                    TScalar const& bottomV,
                    TScalar const& topU,
                    TScalar const& rightBottomV,
                    TScalar const& currentX,
                    TScalar const& currentY,
                    TScalar const& leftX,
                    TScalar const& rightX,
                    TScalar const& bottomY,
                    TScalar const& topY,
                    TScalar const& gamma) {
  return du2dx(currentU,
               leftU,
               rightU,
               currentX,
               leftX,
               rightX,
               gamma)
         + duvdy(currentU,
                 currentV,
                 bottomU,
                 bottomV,
                 rightBottomV,
                 rightV,
                 topU,
                 currentY,
                 currentX,
                 bottomY,
                 topY,
                 rightX,
                 gamma);
}

template <typename TScalar>
inline TScalar
computeConvection3D(TScalar const& currentU,
                    TScalar const& currentV,
                    TScalar const& currentW,
                    TScalar const& leftU,
                    TScalar const& rightU,
                    TScalar const& rightV,
                    TScalar const& rightW,
                    TScalar const& bottomU,
                    TScalar const& bottomV,
                    TScalar const& topU,
                    TScalar const& backU,
                    TScalar const& backW,
                    TScalar const& frontU,
                    TScalar const& rightBottomV,
                    TScalar const& rightBackW,
                    TScalar const& currentX,
                    TScalar const& currentY,
                    TScalar const& currentZ,
                    TScalar const& leftX,
                    TScalar const& rightX,
                    TScalar const& bottomY,
                    TScalar const& topY,
                    TScalar const& backZ,
                    TScalar const& frontZ,
                    TScalar const& gamma) {
  return du2dx(currentU,
               leftU,
               rightU,
               currentX,
               leftX,
               rightX,
               gamma)
         + duvdy(currentU,
                 currentV,
                 bottomU,
                 bottomV,
                 rightBottomV,
                 rightV,
                 topU,
                 currentY,
                 currentX,
                 bottomY,
                 topY,
                 rightX,
                 gamma)
         + duvdy(currentU, // duwdz
                 currentW,
                 backU,
                 backW,
                 rightBackW,
                 rightW,
                 frontU,
                 currentZ,
                 currentX,
                 backZ,
                 frontZ,
                 rightX,
                 gamma);
}

template <int TD>
class DiffusionProcessing {
  template <typename TCellAccessor>
  static inline typename TCellAccessor::VectorDsType
  compute(TCellAccessor const& accessor) {}
};

template <>
class DiffusionProcessing<2> {
public:
  template <typename TCellAccessor>
  static inline
  typename TCellAccessor::VectorDsType
  compute(TCellAccessor const& accessor) {
    typedef typename TCellAccessor::VectorDsType Vector;
    Vector result;

    for (int d1 = 0; d1 < 2; ++d1) {
      int d2 = d1 + 1;

      if (d1 == 1) {
        d2 = 0;
      }

      result(d1) = computeDiffusion2D(
        accessor.velocity(d1),
        accessor.velocity(d1, -1, d1),
        accessor.velocity(d1, +1, d1),
        accessor.velocity(d2, -1, d1),
        accessor.velocity(d2, +1, d1),
        accessor.width(d1),
        accessor.width(d2),
        accessor.width(d1, +1, d1),
        accessor.width(d2, -1, d2),
        accessor.width(d2, +1, d2));
    }

    return result;
  }
};

template <>
class DiffusionProcessing<3> {
public:
  template <typename TCellAccessor>
  static inline
  typename TCellAccessor::VectorDsType
  compute(TCellAccessor const& accessor) {
    typedef typename TCellAccessor::VectorDsType Vector;

    Vector result;

    for (int d1 = 0; d1 < 3; ++d1) {
      int d2 = d1 + 1;
      int d3 = d2 + 1;

      if (d1 == 1) {
        d2 = 0;
        d3 = 2;
      } else if (d1 == 2) {
        d2 = 0;
        d3 = 1;
      }

      result(d1) = computeDiffusion3D(
        accessor.velocity(d1),
        accessor.velocity(d1, -1, d1),
        accessor.velocity(d1, +1, d1),
        accessor.velocity(d2, -1, d1),
        accessor.velocity(d2, +1, d1),
        accessor.velocity(d3, -1, d1),
        accessor.velocity(d3, +1, d1),
        accessor.width(d1),
        accessor.width(d2),
        accessor.width(d3),
        accessor.width(d1, +1, d1),
        accessor.width(d2, -1, d2),
        accessor.width(d2, +1, d2),
        accessor.width(d3, -1, d3),
        accessor.width(d3, +1, d3));
    }

    return result;
  }
};

template <int TD>
class ConvectionProcessing {
  template <typename TCellAccessor,
            typename TSimulationParameters>
  static inline typename TCellAccessor::VectorDsType
  compute(TCellAccessor const&         accessor,
          TSimulationParameters const& simulationParameters) {}
};

template <>
class ConvectionProcessing<2> {
public:
  template <typename TCellAccessor,
            typename TSimulationParameters>
  static inline typename TCellAccessor::VectorDsType
  compute(TCellAccessor const&         accessor,
          TSimulationParameters const* simulationParameters) {
    typedef typename TCellAccessor::VectorDsType Vector;

    Vector result;

    for (int d1 = 0; d1 < 2; ++d1) {
      int d2 = d1 + 1;

      if (d1 == 1) {
        d2 = 0;
      }
      result(d1) = computeConvection2D(
        accessor.velocity(d1),
        accessor.velocity(d2),
        accessor.velocity(d1, -1, d1),
        accessor.velocity(d1, +1, d1),
        accessor.velocity(d1, +1, d2),
        accessor.velocity(d2, -1, d1),
        accessor.velocity(d2, -1, d2),
        accessor.velocity(d2, +1, d1),
        accessor.velocity(d2, -1, d1, +1, d2),
        accessor.width(d1),
        accessor.width(d2),
        accessor.width(d1, -1, d1),
        accessor.width(d1, +1, d1),
        accessor.width(d2, -1, d2),
        accessor.width(d2, +1, d2),
        simulationParameters->gamma());
    }

    return result;
  }
};

template <>
class ConvectionProcessing<3> {
public:
  template <typename TCellAccessor,
            typename TSimulationParameters>
  static inline typename TCellAccessor::VectorDsType
  compute(TCellAccessor const&         accessor,
          TSimulationParameters const* simulationParameters) {
    typedef typename TCellAccessor::VectorDsType Vector;
    Vector result;

    for (int d1 = 0; d1 < 3; ++d1) {
      int d2 = d1 + 1;
      int d3 = d2 + 1;

      if (d1 == 1) {
        d2 = 0;
        d3 = 2;
      } else if (d1 == 2) {
        d2 = 0;
        d3 = 1;
      }

      result(d1) = computeConvection3D(
        accessor.velocity(d1),
        accessor.velocity(d2),
        accessor.velocity(d3),
        accessor.velocity(d1, -1, d1),
        accessor.velocity(d1, +1, d1),
        accessor.velocity(d1, +1, d2),
        accessor.velocity(d1, +1, d3),
        accessor.velocity(d2, -1, d1),
        accessor.velocity(d2, -1, d2),
        accessor.velocity(d2, +1, d1),
        accessor.velocity(d3, -1, d1),
        accessor.velocity(d3, -1, d3),
        accessor.velocity(d3, +1, d1),
        accessor.velocity(d2, -1, d1, +1, d2),
        accessor.velocity(d3, -1, d1, +1, d3),
        accessor.width(d1),
        accessor.width(d2),
        accessor.width(d3),
        accessor.width(d1, -1, d1),
        accessor.width(d1, +1, d1),
        accessor.width(d2, -1, d2),
        accessor.width(d2, +1, d2),
        accessor.width(d3, -1, d3),
        accessor.width(d3, +1, d3),
        simulationParameters->gamma());
    }

    return result;
  }
};

template <int TSolverType>
struct PressureProcessing {};

template <>
struct PressureProcessing<0> {
  template <typename TCellAccessor>
  static inline typename TCellAccessor::VectorDsType
  grad(TCellAccessor const& accessor) {
    typedef typename TCellAccessor::VectorDsType Vector;
    Vector result;

    for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
      result(d)
        = (accessor.pressure(d, +1) - accessor.pressure())
          / (0.5 * (accessor.width(d, +1, d) + accessor.width(d)));
    }

    return result;
  }

  template <typename TCellAccessor, typename TScalar>
  static void
  initialize(TCellAccessor const& accessor, TScalar const& value) {
    accessor.pressure() = value;
  }
};

template <>
struct PressureProcessing<1> {
  template <typename TCellAccessor>
  static inline typename TCellAccessor::VectorDsType
  grad(TCellAccessor const& accessor) {
    typedef typename TCellAccessor::VectorDsType Vector;
    Vector result;

    for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
      result(d)
        = (accessor.projectionTerm(d, +1) - accessor.projectionTerm())
          / (0.5 * (accessor.width(d, +1, d) + accessor.width(d)));
    }

    return result;
  }

  template <typename TCellAccessor, typename TScalar>
  static void
  initialize(TCellAccessor const& accessor, TScalar const& value) {
    accessor.pressure()       = value;
    accessor.projectionTerm() = value;
  }
};

template <typename TCellAccessor,
          typename TParallelDistribution,
          typename TParameters>
class VpeStencilGenerator {
public:
  using CellAccessorType         = TCellAccessor;
  using ParallelDistributionType = TParallelDistribution;
  using ParametersType           = TParameters;
  using ScalarType               = typename CellAccessorType::ScalarType;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

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
  get(TCellAccessor const& accessor,
      TStencil&            stencil,
      TRow&                row,
      TColumns&            columns) const {
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

      auto leftIndex = accessor.relativeIndex(d, -1);
      leftIndex += corner;
      auto rightIndex = accessor.relativeIndex(d, +1);
      rightIndex += corner;

      columns[2 * d].i     = leftIndex(0);
      columns[2 * d].j     = leftIndex(1);
      columns[2 * d + 1].i = rightIndex(0);
      columns[2 * d + 1].j = rightIndex(1);

      if (Dimensions == 3) {
        columns[2 * d].k     = leftIndex(2);
        columns[2 * d + 1].k = rightIndex(2);
      }
    }
    auto currentIndex = accessor.index();
    currentIndex             += corner;
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

      diagonalCoeff            = coeff * diagonalCoeff;
      stencil[2 * Dimensions] += diagonalCoeff;
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

template <typename TCellAccessor,
          int TDimension>
class VpeRhsGenerator {
public:
  using CellAccessorType = TCellAccessor;
  using ScalarType       = typename CellAccessorType::ScalarType;

  inline ScalarType
  get(CellAccessorType const& accessor) const {
    ScalarType result = accessor.fgh(TDimension);

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

template <typename TCellAccessor,
          typename TParallelDistribution>
class PpeStencilGenerator {
public:
  using CellAccessorType         = TCellAccessor;
  using ParallelDistributionType = TParallelDistribution;
  using ScalarType               = typename CellAccessorType::ScalarType;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

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

      auto leftIndex = accessor.relativeIndex(d, -1);
      leftIndex += corner;
      auto rightIndex = accessor.relativeIndex(d, +1);
      rightIndex += corner;

      columns[2 * d].i     = leftIndex(0);
      columns[2 * d].j     = leftIndex(1);
      columns[2 * d + 1].i = rightIndex(0);
      columns[2 * d + 1].j = rightIndex(1);

      if (Dimensions == 3) {
        columns[2 * d].k     = leftIndex(2);
        columns[2 * d + 1].k = rightIndex(2);
      }
    }
    auto currentIndex = accessor.index();
    currentIndex             += corner;
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

template <typename TCellAccessor>
class PpeRhsGenerator {
public:
  using CellAccessorType = TCellAccessor;

  using ScalarType = typename CellAccessorType::ScalarType;

  PpeRhsGenerator() {}

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

    return result;
  }

private:
  ScalarType const* _dt;
};

template <typename TCellAccessor>
class PpeResultAcquirer1 {
public:
  using CellAccessorType = TCellAccessor;

  using ScalarType = typename CellAccessorType::ScalarType;

  inline void
  set(CellAccessorType const& accessor,
      ScalarType const&       value) const {
    accessor.pressure() = value;
  }
};

template <typename TCellAccessor>
class PpeResultAcquirer2 {
public:
  using CellAccessorType = TCellAccessor;
  using ScalarType       = typename CellAccessorType::ScalarType;

  PpeResultAcquirer2() {}

  void
  initialize(ScalarType const* dt) { _dt = dt; }

  inline void
  set(CellAccessorType const& accessor,
      ScalarType const&       value) const {
    // accessor.pressure()       =  accessor.projectionTerm() + value;
    accessor.projectionTerm() = value;
  }

private:
  ScalarType const* _dt;
};
}
}
