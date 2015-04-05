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
             typename GridType::BaseType,
             SfsfdHandlersBuilderTraits::getPressure,
             SfsfdHandlersBuilderTraits::setPressure,
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
