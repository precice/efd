#pragma once

#include "Configuration.hpp"
#include "GhostLayer/InitializationActions.hpp"
#include "GhostLayer/InitializationHandler.hpp"
#include "GhostLayer/MpiExchangeHandler.hpp"
#include "GhostLayer/PetscExchangeActions.hpp"
#include "GhostLayer/PetscExchangeHandler.hpp"
#include "GhostLayer/PressureStencilHanler.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolver,
          int TDimension,
          int TDirection>
class SfsfdGhostHandlersBuilder {
public:
  using SolverType = TSolver;

  enum {
    Dimensions = SolverType::Dimensions
  };

  using CellAccessorType = typename SolverType::CellAccessorType;

  using GridType = typename SolverType::GridType;

  using GhostHandlersType = typename SolverType::GhostHandlersType;

  using ScalarType = typename SolverType::ScalarType;

  using VectorDsType = typename SolverType::VectorDsType;

  template <int TDimension2>
  static ScalarType*
  fghAccessor(CellAccessorType const& accessor) {
    return &accessor.fgh(TDimension2);
  }

  static VectorDsType*
  velocityAccessor(CellAccessorType const& accessor) {
    return &accessor.velocity();
  }

  static ScalarType*
  pressureAccessor(CellAccessorType const& accessor) {
    return &accessor.pressure();
  }

  typedef
    GhostLayer::MpiExchange::Handler
    <ScalarType,
     typename GridType::BaseType,
     SfsfdGhostHandlersBuilder::template fghAccessor<TDimension>,
     TDimension,
     TDirection>
    FghMpiExchangeHandler;

  typedef
    GhostLayer::MpiExchange::Handler
    <ScalarType,
     typename GridType::BaseType,
     SfsfdGhostHandlersBuilder::pressureAccessor,
     TDimension,
     TDirection>
    PressureMpiExchangeHandler;
  typedef
    GhostLayer::MpiExchange::Handler
    <VectorDsType,
     typename GridType::BaseType,
     SfsfdGhostHandlersBuilder::velocityAccessor,
     TDimension,
     TDirection>
    VelocityMpiExchangeHandler;

  typedef
    GhostLayer::Initialization::MovingWallFghAction
    <typename GridType::BaseType,
     ScalarType,
     Dimensions,
     TDimension,
     TDirection>
    MovingWallFghInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::BaseType,
     MovingWallFghInitializationAction,
     Dimensions,
     TDimension,
     TDirection>
    MovingWallFghInitialization;
  typedef
    GhostLayer::Initialization::InputFghAction
    <typename GridType::BaseType,
     ScalarType,
     Dimensions,
     TDimension,
     TDirection>
    InputFghInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::BaseType,
     InputFghInitializationAction,
     Dimensions,
     TDimension,
     TDirection>
    InputFghInitialization;
  typedef
    GhostLayer::Initialization::ParabolicInputFghAction
    <typename GridType::BaseType,
     ScalarType,
     Dimensions,
     TDimension,
     TDirection>
    ParabolicInputFghInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::BaseType,
     ParabolicInputFghInitializationAction,
     Dimensions,
     TDimension,
     TDirection>
    ParabolicInputFghInitialization;
  typedef
    GhostLayer::Initialization::OutputFghAction
    <typename GridType::BaseType,
     ScalarType,
     Dimensions,
     TDimension,
     TDirection>
    OutputFghInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::BaseType,
     OutputFghInitializationAction,
     Dimensions,
     TDimension,
     TDirection>
    OutputFghInitialization;

  typedef
    GhostLayer::LsStencilGenerator::Handler
    <GridType, TDimension, TDirection>
    PpeStencilGenerationHandler;

  typedef
    GhostLayer::PetscExchange::ConstantRhsGenerationAction
    <CellAccessorType, TDimension>
    PpeRhsGenerationAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, PpeRhsGenerationAction, TDimension, TDirection>
    PpeRhsGenerationHandler;

  typedef
    GhostLayer::PetscExchange::PpeRhsAcquiererAction
    PpeRhsAcquiererAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, PpeRhsAcquiererAction, TDimension, TDirection>
    PpeRhsAcquiererHandler;

  typedef
    GhostLayer::Initialization::MovingWallVelocityAction
    <typename GridType::BaseType, ScalarType, Dimensions, TDimension,
     TDirection>
    MovingWallVelocityInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::BaseType, MovingWallVelocityInitializationAction,
     Dimensions,
     TDimension, TDirection>
    MovingWallVelocityInitialization;
  typedef
    GhostLayer::Initialization::InputVelocityAction
    <typename GridType::BaseType, ScalarType, Dimensions, TDimension,
     TDirection>
    InputVelocityInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::BaseType, InputVelocityInitializationAction, Dimensions,
     TDimension, TDirection>
    InputVelocityInitialization;
  typedef
    GhostLayer::Initialization::ParabolicInputVelocityAction
    <typename GridType::BaseType, ScalarType, Dimensions, TDimension,
     TDirection>
    ParabolicInputVelocityInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::BaseType, ParabolicInputVelocityInitializationAction,
     Dimensions,
     TDimension, TDirection>
    ParabolicInputVelocityInitialization;
  typedef
    GhostLayer::Initialization::OutputVelocityAction
    <typename GridType::BaseType, ScalarType, Dimensions, TDimension,
     TDirection>
    OutputVelocityInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::BaseType, OutputVelocityInitializationAction,
     Dimensions,
     TDimension, TDirection>
    OutputVelocityInitialization;

public:
  SfsfdGhostHandlersBuilder(Configuration* configuration,
                            SolverType*    simulation)
    : _configuration(configuration),
    _solver(simulation) {}

  inline void
  setAsInput() {
    if (TDirection == 0) {
      _solver->_ghostHandlers.fghInitialization[TDimension][TDirection] =
        InputFghInitialization::getHandler(
          &_solver->_memory.grid()->indentedBoundaries[TDimension][
            TDirection],
          _solver->_memory.parallelDistribution(),
          new InputFghInitializationAction(_configuration));
    }

    _solver->_ghostHandlers.ppeStencilGeneratorStack[TDimension][TDirection] =
      PpeStencilGenerationHandler::getNeumannMiddle(
        _solver->_memory.grid(),
        _solver->_memory.
                                          parallelDistribution());
    _solver->_ghostHandlers.ppeRhsGeneratorStack[TDimension][TDirection] =
      PpeRhsGenerationHandler::getHandler(_solver->_memory.grid(),
                                          _solver->_memory.
                                          parallelDistribution(),
                                          new PpeRhsGenerationAction(0.0));

    if (TDirection == 1) {
      _solver->_ghostHandlers.ppeRhsAcquiererStack[TDimension][TDirection] =
        PpeRhsAcquiererHandler::getHandler(
          _solver->_memory.grid(),
          _solver->_memory.parallelDistribution(),
          new PpeRhsAcquiererAction());
    }

    _solver->_ghostHandlers.velocityInitialization[TDimension][TDirection] =
      InputVelocityInitialization::getHandler(
        &_solver->_memory.grid()->boundaries[TDimension][TDirection],
        _solver->_memory.parallelDistribution(),
        new InputVelocityInitializationAction(
          _configuration,
          &_solver->_memory.maxVelocity()));
  }

  inline void
  setAsParabolicInput() {
    if (TDirection == 0) {
      _solver->_ghostHandlers.fghInitialization[TDimension][TDirection] =
        ParabolicInputFghInitialization::getHandler(
          &_solver->_memory.grid()->indentedBoundaries[TDimension][
            TDirection],
          _solver->_memory.parallelDistribution(),
          new ParabolicInputFghInitializationAction(_configuration));
    }

    _solver->_ghostHandlers.ppeStencilGeneratorStack[TDimension][TDirection] =
      PpeStencilGenerationHandler::getNeumannMiddle(
        _solver->_memory.grid(),
        _solver->_memory.
                                          parallelDistribution());
    _solver->_ghostHandlers.ppeRhsGeneratorStack[TDimension][TDirection] =
      PpeRhsGenerationHandler::getHandler(_solver->_memory.grid(),
                                          _solver->_memory.
                                          parallelDistribution(),
                                          new PpeRhsGenerationAction(0.0));

    if (TDirection == 1) {
      _solver->_ghostHandlers.ppeRhsAcquiererStack[TDimension][TDirection] =
        PpeRhsAcquiererHandler::getHandler(
          _solver->_memory.grid(),
          _solver->_memory.parallelDistribution(),
          new PpeRhsAcquiererAction());
    }

    _solver->_ghostHandlers.velocityInitialization[TDimension][TDirection] =
      ParabolicInputVelocityInitialization::getHandler(
        &_solver->_memory.grid()->boundaries[TDimension][TDirection],
        _solver->_memory.parallelDistribution(),
        new ParabolicInputVelocityInitializationAction(
          _configuration,
          &_solver->_memory.maxVelocity()));
  }

  inline void
  setAsOutput() {
    if (TDirection == 0) {
      _solver->_ghostHandlers.fghInitialization[TDimension][TDirection] =
        OutputFghInitialization::getHandler(
          &_solver->_memory.grid()->indentedBoundaries[TDimension][
            TDirection],
          _solver->_memory.parallelDistribution(),
          new OutputFghInitializationAction());
    }

    _solver->_ghostHandlers.ppeStencilGeneratorStack[TDimension][TDirection] =
      PpeStencilGenerationHandler::getDirichletMiddle(
        _solver->_memory.grid(),
        _solver->_memory.
                                          parallelDistribution());
    _solver->_ghostHandlers.ppeRhsGeneratorStack[TDimension][TDirection] =
      PpeRhsGenerationHandler::getHandler(_solver->_memory.grid(),
                                          _solver->_memory.
                                          parallelDistribution(),
                                          new PpeRhsGenerationAction(0.0));

    if (TDirection == 1) {
      _solver->_ghostHandlers.ppeRhsAcquiererStack[TDimension][TDirection] =
        PpeRhsAcquiererHandler::getHandler(
          _solver->_memory.grid(),
          _solver->_memory.parallelDistribution(),
          new PpeRhsAcquiererAction());
    }

    _solver->_ghostHandlers.velocityInitialization[TDimension][TDirection]
      = OutputVelocityInitialization::getHandler(
      &_solver->_memory.grid()->boundaries[TDimension][TDirection],
      _solver->_memory.parallelDistribution(),
      new OutputVelocityInitializationAction(
        &_solver->_memory.maxVelocity()));
  }

  void
  setAsMpiExchange() {
    auto leftIndent  = _solver->_memory.grid()->innerGrid.leftIndent();
    auto rightIndent = _solver->_memory.grid()->innerGrid.rightIndent();

    if (TDirection == 0) {
      _solver->_ghostHandlers.mpiFghExchangeStack[TDimension][TDirection] =
        FghMpiExchangeHandler::getReceiveHandler(
          _solver->_memory.grid(),
          _solver->_memory.parallelDistribution(),
          leftIndent,
          rightIndent);
      _solver->_ghostHandlers.mpiPressureExchangeStack[TDimension][TDirection]
        = PressureMpiExchangeHandler::getSendHandler(
        _solver->_memory.grid(),
        _solver->_memory.parallelDistribution(),
        leftIndent,
        rightIndent);

      for (int d = 0; d < Dimensions; ++d) {
        if (d != TDimension) {
          rightIndent(d) = 0;
        }
      }
    } else {
      _solver->_ghostHandlers.mpiFghExchangeStack[TDimension][TDirection] =
        FghMpiExchangeHandler::getSendHandler(
          _solver->_memory.grid(),
          _solver->_memory.parallelDistribution(),
          leftIndent,
          rightIndent);
      _solver->_ghostHandlers.mpiPressureExchangeStack[TDimension][TDirection]
        = PressureMpiExchangeHandler::getReceiveHandler(
        _solver->_memory.grid(),
        _solver->_memory.
        parallelDistribution(),
        leftIndent,
        rightIndent);

      for (int d = 0; d < Dimensions; ++d) {
        if (d != TDimension) {
          leftIndent(d) = 0;
        }
      }
    }

    _solver->_ghostHandlers.mpiVelocityExchangeStack[TDimension][TDirection] =
      VelocityMpiExchangeHandler::getExchangeHandler(
        _solver->_memory.grid(),
        _solver->_memory.parallelDistribution(),
        leftIndent,
        rightIndent);
  }

protected:
  Configuration* _configuration;
  SolverType*    _solver;
};
}
}
