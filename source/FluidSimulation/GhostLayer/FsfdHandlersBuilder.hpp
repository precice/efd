#pragma once

#include "InitializationActions.hpp"
#include "InitializationHandler.hpp"
#include "MpiExchangeHandler.hpp"
#include "PetscExchangeActions.hpp"
#include "PetscExchangeHandler.hpp"
#include "PressureStencilHanler.hpp"

#include "FluidSimulation/Configuration.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <typename THandlersBuilderTraits>
class FsfdHandlersBuilder {
public:
  using SolverTraitsType = typename THandlersBuilderTraits::SolverTraitsType;

  using SolverType = typename SolverTraitsType::SolverType;

  enum {
    Dimension  = THandlersBuilderTraits::Dimension,
    Direction  = THandlersBuilderTraits::Direction,
    Dimensions = SolverTraitsType::Dimensions
  };

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  template <int TDimension2>
  static ScalarType*
  fghAccessor(CellAccessorType const& accessor) {
    return &accessor.fgh(TDimension2);
  }

  static VectorDsType*
  velocityAccessor(CellAccessorType const& accessor) {
    return &accessor.velocity();
  }

  typedef
    MpiExchange::Handler
    <ScalarType,
     typename GridType::BaseType,
     FsfdHandlersBuilder::template fghAccessor<Dimension>,
     Dimension,
     Direction>
    FghMpiExchangeHandler;
  typedef
    typename THandlersBuilderTraits::PressureMpiExchangeHandler
    PressureMpiExchangeHandler;
  typedef
    MpiExchange::Handler
    <VectorDsType,
     typename GridType::BaseType,
     FsfdHandlersBuilder::velocityAccessor,
     Dimension,
     Direction>
    VelocityMpiExchangeHandler;

  typedef
    Initialization::MovingWallFghAction
    <typename GridType::BaseType,
     ScalarType,
     Dimensions,
     Dimension,
     Direction>
    MovingWallFghInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType,
     MovingWallFghInitializationAction,
     Dimensions,
     Dimension,
     Direction>
    MovingWallFghInitialization;
  typedef
    Initialization::InputFghAction
    <typename GridType::BaseType,
     ScalarType,
     Dimensions,
     Dimension,
     Direction>
    InputFghInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType,
     InputFghInitializationAction,
     Dimensions,
     Dimension,
     Direction>
    InputFghInitialization;
  typedef
    Initialization::ParabolicInputFghAction
    <typename GridType::BaseType,
     ScalarType,
     Dimensions,
     Dimension,
     Direction>
    ParabolicInputFghInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType,
     ParabolicInputFghInitializationAction,
     Dimensions,
     Dimension,
     Direction>
    ParabolicInputFghInitialization;
  typedef
    Initialization::OutputFghAction
    <typename GridType::BaseType,
     ScalarType,
     Dimensions,
     Dimension,
     Direction>
    OutputFghInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType,
     OutputFghInitializationAction,
     Dimensions,
     Dimension,
     Direction>
    OutputFghInitialization;

  typedef
    LsStencilGenerator::Handler
    <GridType, Dimension, Direction>
    PpeStencilGenerationHandler;

  typedef
    PetscExchange::ConstantRhsGenerationAction
    <CellAccessorType, Dimension>
    PpeRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, PpeRhsGenerationAction, Dimension, Direction>
    PpeRhsGenerationHandler;

  typedef
    typename THandlersBuilderTraits::PpeRhsAcquiererAction
    PpeRhsAcquiererAction;
  typedef
    PetscExchange::Handler
    <GridType, PpeRhsAcquiererAction, Dimension, Direction>
    PpeRhsAcquiererHandler;

  typedef
    Initialization::MovingWallVelocityAction
    <typename GridType::BaseType, ScalarType, Dimensions, Dimension,
     Direction>
    MovingWallVelocityInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType, MovingWallVelocityInitializationAction,
     Dimensions,
     Dimension, Direction>
    MovingWallVelocityInitialization;
  typedef
    Initialization::InputVelocityAction
    <typename GridType::BaseType, ScalarType, Dimensions, Dimension,
     Direction>
    InputVelocityInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType, InputVelocityInitializationAction, Dimensions,
     Dimension, Direction>
    InputVelocityInitialization;
  typedef
    Initialization::ParabolicInputVelocityAction
    <typename GridType::BaseType, ScalarType, Dimensions, Dimension,
     Direction>
    ParabolicInputVelocityInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType, ParabolicInputVelocityInitializationAction,
     Dimensions,
     Dimension, Direction>
    ParabolicInputVelocityInitialization;
  typedef
    Initialization::OutputVelocityAction
    <typename GridType::BaseType, ScalarType, Dimensions, Dimension,
     Direction>
    OutputVelocityInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType, OutputVelocityInitializationAction,
     Dimensions,
     Dimension, Direction>
    OutputVelocityInitialization;

public:
  FsfdHandlersBuilder(Configuration* configuration,
                      SolverType*    simulation)
    : _configuration(configuration),
    _solver(simulation) {}

  inline void
  setAsInput() {
    if (Direction == 0) {
      _solver->ghostHandlers()->fghInitialization[Dimension][Direction]
        = InputFghInitialization::getHandler(
        &_solver->memory()->grid()->indentedBoundaries[Dimension][
          Direction],
        _solver->memory()->parallelDistribution(),
        new InputFghInitializationAction(_configuration));
    }

    _solver->ghostHandlers()->ppeStencilGeneratorStack[Dimension][Direction]
      = PpeStencilGenerationHandler::getNeumannMiddle(
      _solver->memory()->grid(),
      _solver->memory()->parallelDistribution());
    _solver->ghostHandlers()->ppeRhsGeneratorStack[Dimension][Direction]
      = PpeRhsGenerationHandler::getHandler(_solver->memory()->grid(),
                                            _solver->memory()->
                                            parallelDistribution(),
                                            new PpeRhsGenerationAction(0.0));

    if (Direction == 1) {
      _solver->ghostHandlers()->ppeRhsAcquiererStack[Dimension][Direction]
        = PpeRhsAcquiererHandler::getHandler(
        _solver->memory()->grid(),
        _solver->memory()->parallelDistribution(),
        new PpeRhsAcquiererAction());
    }

    _solver->ghostHandlers()->velocityInitialization[Dimension][Direction]
      = InputVelocityInitialization::getHandler(
      &_solver->memory()->grid()->boundaries[Dimension][Direction],
      _solver->memory()->parallelDistribution(),
      new InputVelocityInitializationAction(
        _configuration,
        &_solver->memory()->maxVelocity()));
  }

  inline void
  setAsParabolicInput() {
    if (Direction == 0) {
      _solver->ghostHandlers()->fghInitialization[Dimension][Direction]
        = ParabolicInputFghInitialization::getHandler(
        &_solver->memory()->grid()
        ->indentedBoundaries[Dimension][Direction],
        _solver->memory()->parallelDistribution(),
        new ParabolicInputFghInitializationAction(_configuration));
    }

    _solver->ghostHandlers()->ppeStencilGeneratorStack[Dimension][Direction]
      = PpeStencilGenerationHandler::getNeumannMiddle(
      _solver->memory()->grid(),
      _solver->memory()->parallelDistribution());
    _solver->ghostHandlers()->ppeRhsGeneratorStack[Dimension][Direction]
      = PpeRhsGenerationHandler::getHandler(
      _solver->memory()->grid(),
      _solver->memory()->parallelDistribution(),
      new PpeRhsGenerationAction(0.0));

    if (Direction == 1) {
      _solver->ghostHandlers()->ppeRhsAcquiererStack[Dimension][Direction]
        = PpeRhsAcquiererHandler::getHandler(
        _solver->memory()->grid(),
        _solver->memory()->parallelDistribution(),
        new PpeRhsAcquiererAction());
    }

    _solver->ghostHandlers()->velocityInitialization[Dimension][Direction]
      = ParabolicInputVelocityInitialization::getHandler(
      &_solver->memory()->grid()->boundaries[Dimension][Direction],
      _solver->memory()->parallelDistribution(),
      new ParabolicInputVelocityInitializationAction(
        _configuration,
        &_solver->memory()->maxVelocity()));
  }

  inline void
  setAsOutput() {
    if (Direction == 0) {
      _solver->ghostHandlers()->fghInitialization[Dimension][Direction]
        = OutputFghInitialization::getHandler(
        &_solver->memory()->grid()
        ->indentedBoundaries[Dimension][Direction],
        _solver->memory()->parallelDistribution(),
        new OutputFghInitializationAction());
    }

    _solver->ghostHandlers()->ppeStencilGeneratorStack[Dimension][Direction]
      = PpeStencilGenerationHandler::getDirichletMiddle(
      _solver->memory()->grid(),
      _solver->memory()->parallelDistribution());
    _solver->ghostHandlers()->ppeRhsGeneratorStack[Dimension][Direction]
      = PpeRhsGenerationHandler::getHandler(
      _solver->memory()->grid(),
      _solver->memory()->parallelDistribution(),
      new PpeRhsGenerationAction(0.0));

    if (Direction == 1) {
      _solver->ghostHandlers()->ppeRhsAcquiererStack[Dimension][Direction]
        = PpeRhsAcquiererHandler::getHandler(
        _solver->memory()->grid(),
        _solver->memory()->parallelDistribution(),
        new PpeRhsAcquiererAction());
    }

    _solver->ghostHandlers()->velocityInitialization[Dimension][Direction]
      = OutputVelocityInitialization::getHandler(
      &_solver->memory()->grid()->boundaries[Dimension][Direction],
      _solver->memory()->parallelDistribution(),
      new OutputVelocityInitializationAction(
        &_solver->memory()->maxVelocity()));
  }

  void
  setAsMpiExchange() {
    auto leftIndent  = _solver->memory()->grid()->innerGrid.leftIndent();
    auto rightIndent = _solver->memory()->grid()->innerGrid.rightIndent();

    if (Direction == 0) {
      _solver->ghostHandlers()->mpiFghExchangeStack[Dimension][Direction]
        = FghMpiExchangeHandler::getReceiveHandler(
        _solver->memory()->grid(),
        _solver->memory()->parallelDistribution(),
        leftIndent,
        rightIndent);

      _solver->ghostHandlers()->mpiPressureExchangeStack[Dimension][Direction]
        = PressureMpiExchangeHandler::getSendHandler(
        _solver->memory()->grid(),
        _solver->memory()->parallelDistribution(),
        leftIndent,
        rightIndent);

      for (int d = 0; d < Dimensions; ++d) {
        if (d != Dimension) {
          rightIndent(d) = 0;
        }
      }
    } else {
      _solver->ghostHandlers()->mpiFghExchangeStack[Dimension][Direction]
        = FghMpiExchangeHandler::getSendHandler(
        _solver->memory()->grid(),
        _solver->memory()->parallelDistribution(),
        leftIndent,
        rightIndent);

      _solver->ghostHandlers()->mpiPressureExchangeStack[Dimension][Direction]
        = PressureMpiExchangeHandler::getReceiveHandler(
        _solver->memory()->grid(),
        _solver->memory()->parallelDistribution(),
        leftIndent,
        rightIndent);

      for (int d = 0; d < Dimensions; ++d) {
        if (d != Dimension) {
          leftIndent(d) = 0;
        }
      }
    }

    _solver->ghostHandlers()->mpiVelocityExchangeStack[Dimension][Direction]
      = VelocityMpiExchangeHandler::getExchangeHandler(
      _solver->memory()->grid(),
      _solver->memory()->parallelDistribution(),
      leftIndent,
      rightIndent);
  }

protected:
  Configuration* _configuration;
  SolverType*    _solver;
};
}
}
}
