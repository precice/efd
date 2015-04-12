#ifndef Fluid_Simulation_GhostLayer_SfsfdHandlersBuilder_hpp
#define Fluid_Simulation_GhostLayer_SfsfdHandlersBuilder_hpp

#include "FsfdHandlersBuilder.hpp"

#include "Simulation/FsfdSolver.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace PetscExchange {
class PpeRhsAcquiererAction1;
}
template <typename TSolverTraits,
          int TDimension,
          int TDirection>
struct SfsfdHandlersBuilderTraits {
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimension  = TDimension,
    Direction  = TDirection,
    Dimensions = SolverTraitsType::Dimensions
  };

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using BaseGridType = typename SolverTraitsType::BaseGridType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  static ScalarType*
  getPressure(CellAccessorType const& accessor) {
    return &accessor.pressure();
  }

  static void
  setPressure(CellAccessorType const& accessor,
              int const&              index,
              ScalarType const&       value) {
    ((void)index);
    accessor.pressure() = value;
  }

  using PressureMpiExchangeHandler
          = MpiExchange::Handler
            <ScalarType,
             1,
             BaseGridType,
             SfsfdHandlersBuilderTraits::getPressure,
             SfsfdHandlersBuilderTraits::setPressure,
             TDimension,
             TDirection>;

  using PpeRhsAcquiererAction
          = PetscExchange::PpeRhsAcquiererAction1;
};

extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 1 >>;

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class SfsfdHandlersBuilder :
  public FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits < TSolverTraits, TDimension, TDirection >> {
public:
  using HandlersBuilderTraitsType
          = SfsfdHandlersBuilderTraits<TSolverTraits, TDimension, TDirection>;

  using BaseType = FsfdHandlersBuilder<HandlersBuilderTraitsType>;

  using SolverTraitsType = TSolverTraits;

  using SolverType = typename SolverTraitsType::SolverType;

  enum {
    Dimension  = TDimension,
    Direction  = TDirection,
    Dimensions = SolverTraitsType::Dimensions
  };

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

public:
  SfsfdHandlersBuilder(Configuration*, SolverType*);

  void
  setAsInput();

  void
  setAsParabolicInput();

  void
  setAsOutput();

  void
  setAsMpiExchange();

protected:
  Configuration* _configuration;
  SolverType*    _solver;
};

extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 1>;

extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 1>;

extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 1>;

extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 1>;

extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 1>;

extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 1>;

extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 1>;

extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 1>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 0>;
extern template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 1>;
}
}
}
#endif
