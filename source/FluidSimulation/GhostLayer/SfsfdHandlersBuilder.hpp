#pragma once

#include "FsfdHandlersBuilder.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
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

  using GridType = typename SolverTraitsType::GridType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  static ScalarType*
  pressureAccessor(CellAccessorType const& accessor) {
    return &accessor.pressure();
  }

  using PressureMpiExchangeHandler
          = MpiExchange::Handler
            <ScalarType,
             1,
             typename GridType::BaseType,
             SfsfdHandlersBuilderTraits::pressureAccessor,
             TDimension,
             TDirection>;

  using PpeRhsAcquiererAction
          = PetscExchange::PpeRhsAcquiererAction1;
};

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

  using PressureMpiExchangeHandler
          = typename HandlersBuilderTraitsType::PressureMpiExchangeHandler;

  using PpeRhsAcquiererAction
          = typename HandlersBuilderTraitsType::PpeRhsAcquiererAction;

public:
  SfsfdHandlersBuilder(Configuration* configuration,
                       SolverType*    simulation)
    : BaseType(configuration, simulation) {}

  inline void
  setAsInput() {
    this->BaseType::setAsInput();
  }

  inline void
  setAsParabolicInput() {
    this->BaseType::setAsParabolicInput();
  }

  inline void
  setAsOutput() {
    this->BaseType::setAsOutput();
  }

  void
  setAsMpiExchange() {
    this->BaseType::setAsMpiExchange();
  }

protected:
  Configuration* _configuration;
  SolverType*    _solver;
};
}
}
}
