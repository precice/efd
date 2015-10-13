#ifndef Fluid_Simulation_GhostLayer_IfsfdHandlersBuilder_hpp
#define Fluid_Simulation_GhostLayer_IfsfdHandlersBuilder_hpp

#include "FsfdHandlersBuilder.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace PetscExchange {
class PpeRhsAcquiererAction2;
}
template <typename TSolverTraits,
          int TDimension,
          int TDirection>
struct IfsfdHandlersBuilderTraits {
  using SolverTraitsType = TSolverTraits;
  enum {
    Dimension  = TDimension,
    Direction  = TDirection,
    Dimensions = SolverTraitsType::Dimensions
  };

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using SubgridType = typename SolverTraitsType::SubgridType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  inline static ScalarType*
  getProjectionTerm(CellAccessorType const& accessor) {
    return &accessor.projectionTerm();
  }

  inline static void
  setProjectionTerm(CellAccessorType const& accessor,
                    int const&              index,
                    ScalarType const&       value) {
    ((void)index);
    accessor.pressure()      += value;
    accessor.projectionTerm() = value;
  }

  using PressureMpiExchangeHandler
          = MpiExchange::Handler
            <ScalarType,
             1,
             SubgridType,
             IfsfdHandlersBuilderTraits::getProjectionTerm,
             IfsfdHandlersBuilderTraits::setProjectionTerm,
             TDimension,
             TDirection>;

  using PpeRhsAcquiererAction
          = PetscExchange::PpeRhsAcquiererAction2;
};

extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 2>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 2>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 2>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 2>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 2>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 2>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 2>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 2>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 2>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 2>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 2>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 2>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 3>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 3>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 3>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 3>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 3>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<0, double, 3>, 2, 1 >>;

extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 3>, 0, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 3>, 0, 1 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 3>, 1, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 3>, 1, 1 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 3>, 2, 0 >>;
extern template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<1, double, 3>, 2, 1 >>;

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class IfsfdHandlersBuilder :
  public FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits < TSolverTraits, TDimension, TDirection >> {
public:
  using HandlersBuilderTraitsType
          = IfsfdHandlersBuilderTraits<TSolverTraits, TDimension, TDirection>;

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
  IfsfdHandlersBuilder(Configuration*, SolverType*);

  void
  setAsInput();

  void
  setAsParabolicInput();

  void
  setAsOutput();

  void
  setAsMpiExchange();
};

extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 0, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 0, 1>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 1, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 1, 1>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 2, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 2, 1>;

extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 0, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 0, 1>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 1, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 1, 1>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 2, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 2, 1>;

extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 0, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 0, 1>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 1, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 1, 1>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 2, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 2, 1>;

extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 0, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 0, 1>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 1, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 1, 1>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 2, 0>;
extern template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 2, 1>;
}
}
}
#endif
