#pragma once

#include "SfsfdGhostHandlersBuilder.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolver,
          int TDimension,
          int TDirection>
class IfsfdGhostHandlersBuilder
: public SfsfdGhostHandlersBuilder<TSolver, TDimension, TDirection>{
public:
 using BaseType
     = SfsfdGhostHandlersBuilder<TSolver, TDimension, TDirection>;

  using SolverType = typename BaseType::SolverType;

  enum {
    Dimensions = BaseType::Dimensions
  };

  using CellAccessorType = typename BaseType::CellAccessorType;

  using GridType = typename BaseType::GridType;

  using GhostHandlersType = typename BaseType::GhostHandlersType;

  using VectorDsType = typename BaseType::VectorDsType;

  typedef
    GhostLayer::MpiExchange::Handler
    <ScalarType,
     typename GridType::Base,
     template fghAccessor<TDimension>,
     TDimension,
     TDirection>
    FghMpiExchangeHandler;
  typedef
    GhostLayer::MpiExchange::Handler
    <ScalarType,
     typename GridType::Base,
     pressureAccessor,
     TDimension,
     TDirection>
    PressureMpiExchangeHandler;
  typedef
    GhostLayer::MpiExchange::Handler
    <VectorDsType,
     typename GridType::Base,
     velocityAccessor,
     TDimension,
     TDirection>
    VelocityMpiExchangeHandler;

  typedef
    GhostLayer::Initialization::MovingWallFghAction
    <typename GridType::Base,
     ScalarType,
     Dimensions,
     TDimension,
     TDirection>
    MovingWallFghInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::Base,
     MovingWallFghInitializationAction,
     Dimensions,
     TDimension,
     TDirection>
    MovingWallFghInitialization;
  typedef
    GhostLayer::Initialization::InputFghAction
    <typename GridType::Base,
     ScalarType,
     Dimensions,
     TDimension,
     TDirection>
    InputFghInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::Base,
     InputFghInitializationAction,
     Dimensions,
     TDimension,
     TDirection>
    InputFghInitialization;
  typedef
    GhostLayer::Initialization::ParabolicInputFghAction
    <typename GridType::Base,
     ScalarType,
     Dimensions,
     TDimension,
     TDirection>
    ParabolicInputFghInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::Base,
     ParabolicInputFghInitializationAction,
     Dimensions,
     TDimension,
     TDirection>
    ParabolicInputFghInitialization;
  typedef
    GhostLayer::Initialization::OutputFghAction
    <typename GridType::Base,
     ScalarType,
     Dimensions,
     TDimension,
     TDirection>
    OutputFghInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::Base,
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
    GhostLayer::LsStencilGenerator::Handler
    <GridType, TDimension, TDirection>
    VpeStencilGenerationHandler;

  typedef
    GhostLayer::PetscExchange::ConstantRhsGenerationAction
    <CellAccessorType, TDimension>
    VpeConstantRhsGenerationAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VpeConstantRhsGenerationAction, TDimension, TDirection>
    VpeConstantRhsGenerationHandler;

  typedef
    GhostLayer::PetscExchange::VpeInputRhsGenerationAction
    <0, TDimension, TDirection>
    VxpeInputRhsGenerationAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VxpeInputRhsGenerationAction, TDimension, TDirection>
    VxpeInputRhsGenerationHandler;

  typedef
    GhostLayer::PetscExchange::
    VpeParabolicInputRhsGenerationAction
    <0, TDimension, TDirection>
    VxpeParabolicInputRhsGenerationAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VxpeParabolicInputRhsGenerationAction, TDimension, TDirection>
    VxpeParabolicInputRhsGenerationHandler;

  typedef
    GhostLayer::PetscExchange::VpeRhsAcquiererAction<0>
    VxpeRhsAcquiererAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VxpeRhsAcquiererAction, TDimension, TDirection>
    VxpeRhsAcquiererHandler;

  typedef
    GhostLayer::PetscExchange::VpeInputRhsGenerationAction
    <1, TDimension, TDirection>
    VypeInputRhsGenerationAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VypeInputRhsGenerationAction, TDimension, TDirection>
    VypeInputRhsGenerationHandler;

  typedef
    GhostLayer::PetscExchange::
    VpeParabolicInputRhsGenerationAction
    <1, TDimension, TDirection>
    VypeParabolicInputRhsGenerationAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VypeParabolicInputRhsGenerationAction, TDimension, TDirection>
    VypeParabolicInputRhsGenerationHandler;

  typedef
    GhostLayer::PetscExchange::VpeRhsAcquiererAction<0>
    VypeRhsAcquiererAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VypeRhsAcquiererAction, TDimension, TDirection>
    VypeRhsAcquiererHandler;

  typedef
    GhostLayer::PetscExchange::VpeInputRhsGenerationAction
    <2, TDimension, TDirection>
    VzpeInputRhsGenerationAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VzpeInputRhsGenerationAction, TDimension, TDirection>
    VzpeInputRhsGenerationHandler;

  typedef
    GhostLayer::PetscExchange::
    VpeParabolicInputRhsGenerationAction
    <2, TDimension, TDirection>
    VzpeParabolicInputRhsGenerationAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VzpeParabolicInputRhsGenerationAction, TDimension, TDirection>
    VzpeParabolicInputRhsGenerationHandler;

  typedef
    GhostLayer::PetscExchange::VpeRhsAcquiererAction<0>
    VzpeRhsAcquiererAction;
  typedef
    GhostLayer::PetscExchange::Handler
    <GridType, VzpeRhsAcquiererAction, TDimension, TDirection>
    VzpeRhsAcquiererHandler;

  typedef
    GhostLayer::Initialization::MovingWallVelocityAction
    <typename GridType::Base, ScalarType, Dimensions, TDimension, TDirection>
    MovingWallVelocityInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::Base, MovingWallVelocityInitializationAction,
     Dimensions,
     TDimension, TDirection>
    MovingWallVelocityInitialization;
  typedef
    GhostLayer::Initialization::InputVelocityAction
    <typename GridType::Base, ScalarType, Dimensions, TDimension, TDirection>
    InputVelocityInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::Base, InputVelocityInitializationAction, Dimensions,
     TDimension, TDirection>
    InputVelocityInitialization;
  typedef
    GhostLayer::Initialization::ParabolicInputVelocityAction
    <typename GridType::Base, ScalarType, Dimensions, TDimension, TDirection>
    ParabolicInputVelocityInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::Base, ParabolicInputVelocityInitializationAction,
     Dimensions,
     TDimension, TDirection>
    ParabolicInputVelocityInitialization;
  typedef
    GhostLayer::Initialization::OutputVelocityAction
    <typename GridType::Base, ScalarType, Dimensions, TDimension, TDirection>
    OutputVelocityInitializationAction;
  typedef
    GhostLayer::Initialization::Handler
    <typename GridType::Base, OutputVelocityInitializationAction, Dimensions,
     TDimension, TDirection>
    OutputVelocityInitialization;

public:
  GhostHandlers(Configuration * configuration,
                SimulationType *                 simulation)
    : _configuration(configuration),
      _simulation(simulation) {}

  inline void
  setAsInput() {
    if (TDirection == 0) {
      _simulation->_handlers->fghInitialization[TDimension][TDirection] =
        InputFghInitialization::getHandler(
          &_simulation->_memory.grid()->indentedBoundaries[TDimension][
            TDirection],
          _simulation->_memory.parallelDistribution(),
          new InputFghInitializationAction(_configuration));
    }

    _simulation->_handlers->ppeStencilGeneratorStack[TDimension][TDirection] =
      PpeStencilGenerationHandler::getNeumannMiddle(
        _simulation->_memory.grid(),
        _simulation->_memory.
                                          parallelDistribution());
    _simulation->_handlers->ppeRhsGeneratorStack[TDimension][TDirection] =
      PpeRhsGenerationHandler::getHandler(_simulation->_memory.grid(),
                                          _simulation->_memory.
                                          parallelDistribution(),
                                          new PpeRhsGenerationAction(0.0));

    if (TDirection == 1) {
      _simulation->_handlers->ppeRhsAcquiererStack[TDimension][TDirection] =
        PpeRhsAcquiererHandler::getHandler(
          _simulation->_memory.grid(),
          _simulation->_memory.parallelDistribution(),
          new PpeRhsAcquiererAction());
    }

    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          _simulation->_handlers->vpeStencilGeneratorStack[d][TDimension][
            TDirection] =
            VpeStencilGenerationHandler::getDirichletLeft(
              _simulation->_memory.grid(),
              _simulation->_memory.parallelDistribution());
        } else {
          _simulation->_handlers->vpeStencilGeneratorStack[d][TDimension][
            TDirection] =
            VpeStencilGenerationHandler::getDirichletRight(
              _simulation->_memory.grid(),
              _simulation->_memory.parallelDistribution());
        }
      } else {
        _simulation->_handlers->vpeStencilGeneratorStack[d][TDimension][
          TDirection] =
          VpeStencilGenerationHandler::getDirichletMiddle(
            _simulation->_memory.grid(),
            _simulation->_memory.parallelDistribution());
      }
    }

    _simulation->_handlers->vpeRhsGeneratorStack[0][TDimension][TDirection] =
      VxpeInputRhsGenerationHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VxpeInputRhsGenerationAction(_configuration));

    if (TDimension == 0 && TDirection == 0) {
      _simulation->_handlers->vpeRhsAcquiererStack[0][TDimension][TDirection]
        =
          VxpeRhsAcquiererHandler::getHandler(
            _simulation->_memory.grid(),
            _simulation->_memory.parallelDistribution(),
            new VxpeRhsAcquiererAction());
    }

    _simulation->_handlers->vpeRhsGeneratorStack[1][TDimension][TDirection] =
      VypeInputRhsGenerationHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VypeInputRhsGenerationAction(_configuration));

    if (TDimension == 1 && TDirection == 0) {
      _simulation->_handlers->vpeRhsAcquiererStack[1][TDimension][TDirection]
        =
          VypeRhsAcquiererHandler::getHandler(
            _simulation->_memory.grid(),
            _simulation->_memory.parallelDistribution(),
            new VypeRhsAcquiererAction());
    }

    if (Dimensions > 2) {
      _simulation->_handlers->vpeRhsGeneratorStack[2][TDimension][TDirection]
        =
          VzpeInputRhsGenerationHandler::getHandler(
            _simulation->_memory.grid(),
            _simulation->_memory.parallelDistribution(),
            new VzpeInputRhsGenerationAction(_configuration));

      if (TDimension == 2 && TDirection == 0) {
        _simulation->_handlers->vpeRhsAcquiererStack[2][TDimension][TDirection
        ] =
          VzpeRhsAcquiererHandler::getHandler(
            _simulation->_memory.grid(),
            _simulation->_memory.parallelDistribution(),
            new VzpeRhsAcquiererAction());
      }
    }

    _simulation->_handlers->velocityInitialization[TDimension][TDirection] =
      InputVelocityInitialization::getHandler(
        &_simulation->_memory.grid()->boundaries[TDimension][TDirection],
        _simulation->_memory.parallelDistribution(),
        new InputVelocityInitializationAction(_configuration, _maxVelocity));
  }

  inline void
  setAsParabolicInput() {
    if (TDirection == 0) {
      _simulation->_handlers->fghInitialization[TDimension][TDirection] =
        ParabolicInputFghInitialization::getHandler(
          &_simulation->_memory.grid()->indentedBoundaries[TDimension][
            TDirection],
          _simulation->_memory.parallelDistribution(),
          new ParabolicInputFghInitializationAction(_configuration));
    }

    _simulation->_handlers->ppeStencilGeneratorStack[TDimension][TDirection] =
      PpeStencilGenerationHandler::getNeumannMiddle(
        _simulation->_memory.grid(),
        _simulation->_memory.
                                          parallelDistribution());
    _simulation->_handlers->ppeRhsGeneratorStack[TDimension][TDirection] =
      PpeRhsGenerationHandler::getHandler(_simulation->_memory.grid(),
                                          _simulation->_memory.
                                          parallelDistribution(),
                                          new PpeRhsGenerationAction(0.0));

    if (TDirection == 1) {
      _simulation->_handlers->ppeRhsAcquiererStack[TDimension][TDirection] =
        PpeRhsAcquiererHandler::getHandler(
          _simulation->_memory.grid(),
          _simulation->_memory.parallelDistribution(),
          new PpeRhsAcquiererAction());
    }

    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          _simulation->_handlers->vpeStencilGeneratorStack[d][TDimension][
            TDirection]
            = VpeStencilGenerationHandler::getDirichletLeft(
            _simulation->_memory.grid(),
            _simulation->_memory.parallelDistribution());
        } else {
          _simulation->_handlers->vpeStencilGeneratorStack[d][TDimension][
            TDirection]
            = VpeStencilGenerationHandler::getDirichletRight(
            _simulation->_memory.grid(),
            _simulation->
            _memory.parallelDistribution());
        }
      } else {
        _simulation->_handlers->vpeStencilGeneratorStack[d][TDimension][
          TDirection]
          = VpeStencilGenerationHandler::getDirichletMiddle(
          _simulation->_memory.grid(),
          _simulation->_memory.parallelDistribution());
      }
    }

    _simulation->_handlers->vpeRhsGeneratorStack[0][TDimension][TDirection] =
      VxpeParabolicInputRhsGenerationHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VxpeParabolicInputRhsGenerationAction(_configuration));

    if (TDimension == 0 && TDirection == 0) {
      _simulation->_handlers->vpeRhsAcquiererStack[0][TDimension][TDirection]
        =
          VxpeRhsAcquiererHandler::getHandler(
            _simulation->_memory.grid(),
            _simulation->_memory.parallelDistribution(),
            new VxpeRhsAcquiererAction());
    }

    _simulation->_handlers->vpeRhsGeneratorStack[1][TDimension][TDirection] =
      VypeParabolicInputRhsGenerationHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VypeParabolicInputRhsGenerationAction(_configuration));

    if (TDimension == 1 && TDirection == 0) {
      _simulation->_handlers->vpeRhsAcquiererStack[1][TDimension][TDirection]
        = VypeRhsAcquiererHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VypeRhsAcquiererAction());
    }

    if (Dimensions > 2) {
      _simulation->_handlers->vpeRhsGeneratorStack[2][TDimension][TDirection]
        = VzpeParabolicInputRhsGenerationHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VzpeParabolicInputRhsGenerationAction(_configuration));

      if (TDimension == 2 && TDirection == 0) {
        _simulation->_handlers
        ->vpeRhsAcquiererStack[2][TDimension][TDirection]
          = VzpeRhsAcquiererHandler::getHandler(
          _simulation->_memory.grid(),
          _simulation->_memory.parallelDistribution(),
          new VzpeRhsAcquiererAction());
      }
    }

    _simulation->_handlers->velocityInitialization[TDimension][TDirection] =
      ParabolicInputVelocityInitialization::getHandler(
        &_simulation->_memory.grid()->boundaries[TDimension][TDirection],
        _simulation->_memory.parallelDistribution(),
        new ParabolicInputVelocityInitializationAction(_configuration,
                                                       _maxVelocity));
  }

  inline void
  setAsOutput() {
    if (TDirection == 0) {
      _simulation->_handlers->fghInitialization[TDimension][TDirection] =
        OutputFghInitialization::getHandler(
          &_simulation->_memory.grid()->indentedBoundaries[TDimension][
            TDirection],
          _simulation->_memory.parallelDistribution(),
          new OutputFghInitializationAction());
    }

    _simulation->_handlers->ppeStencilGeneratorStack[TDimension][TDirection] =
      PpeStencilGenerationHandler::getDirichletMiddle(
        _simulation->_memory.grid(),
        _simulation->_memory.
                                          parallelDistribution());
    _simulation->_handlers->ppeRhsGeneratorStack[TDimension][TDirection] =
      PpeRhsGenerationHandler::getHandler(_simulation->_memory.grid(),
                                          _simulation->_memory.
                                          parallelDistribution(),
                                          new PpeRhsGenerationAction(0.0));

    if (TDirection == 1) {
      _simulation->_handlers->ppeRhsAcquiererStack[TDimension][TDirection] =
        PpeRhsAcquiererHandler::getHandler(
          _simulation->_memory.grid(),
          _simulation->_memory.parallelDistribution(),
          new PpeRhsAcquiererAction());
    }

    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          _simulation->_handlers->vpeStencilGeneratorStack[d][TDimension][
            TDirection]
            = VpeStencilGenerationHandler::getNeumannLeft(
            _simulation->_memory.grid(),
            _simulation->_memory.parallelDistribution());
        } else {
          _simulation->_handlers->vpeStencilGeneratorStack[d][TDimension][
            TDirection]
            = VpeStencilGenerationHandler::getNeumannRight(
            _simulation->_memory.grid(),
            _simulation->_memory.parallelDistribution());
        }
      } else {
        _simulation->_handlers->vpeStencilGeneratorStack[d][TDimension][
          TDirection]
          = VpeStencilGenerationHandler::getNeumannMiddle(
          _simulation->_memory.grid(),
          _simulation->_memory.parallelDistribution());
      }
    }

    int offset = 0;

    if (TDimension == 0) {
      if (TDirection == 1) {
        offset = 1;
      }
    }
    _simulation->_handlers->vpeRhsGeneratorStack[0][TDimension][TDirection] =
      VpeConstantRhsGenerationHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VpeConstantRhsGenerationAction(0.0, offset));

    if (TDimension == 0 && TDirection == 0) {
      _simulation->_handlers->vpeRhsAcquiererStack[0][TDimension][TDirection]
        = VxpeRhsAcquiererHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VxpeRhsAcquiererAction());
    }

    offset = 0;

    if (TDimension == 1) {
      if (TDirection == 1) {
        offset = 1;
      }
    }

    _simulation->_handlers->vpeRhsGeneratorStack[1][TDimension][TDirection] =
      VpeConstantRhsGenerationHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VpeConstantRhsGenerationAction(0.0, offset));

    if (TDimension == 1 && TDirection == 0) {
      _simulation->_handlers->vpeRhsAcquiererStack[1][TDimension][TDirection]
        = VypeRhsAcquiererHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VypeRhsAcquiererAction());
    }

    if (Dimensions > 2) {
      offset = 0;

      if (TDimension == 2) {
        if (TDirection == 1) {
          offset = 1;
        }
      }
      _simulation->_handlers->vpeRhsGeneratorStack[2][TDimension][TDirection]
        = VpeConstantRhsGenerationHandler::getHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        new VpeConstantRhsGenerationAction(0.0, offset));

      if (TDimension == 2 && TDirection == 0) {
        _simulation->_handlers
        ->vpeRhsAcquiererStack[2][TDimension][TDirection]
          = VzpeRhsAcquiererHandler::getHandler(
          _simulation->_memory.grid(),
          _simulation->_memory.parallelDistribution(),
          new VzpeRhsAcquiererAction());
      }
    }

    _simulation->_handlers->velocityInitialization[TDimension][TDirection]
      = OutputVelocityInitialization::getHandler(
      &_simulation->_memory.grid()->boundaries[TDimension][TDirection],
      _simulation->_memory.parallelDistribution(),
      new OutputVelocityInitializationAction(_maxVelocity));
  }

  void
  setAsMpiExchange() {
    auto leftIndent  = _simulation->_memory.grid()->innerGrid.leftIndent();
    auto rightIndent = _simulation->_memory.grid()->innerGrid.rightIndent();

    if (TDirection == 0) {
      _simulation->_handlers->mpiFghExchangeStack[TDimension][TDirection] =
        FghMpiExchangeHandler::getReceiveHandler(_simulation->_memory.grid(),
                                                 _simulation->_memory.
                                                 parallelDistribution(),
                                                 leftIndent,
                                                 rightIndent);
      _simulation->_handlers->mpiPressureExchangeStack[TDimension][TDirection]
        = PressureMpiExchangeHandler::getSendHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        leftIndent,
        rightIndent);

      for (int d = 0; d < Dimensions; ++d) {
        if (d != TDimension) {
          rightIndent(d) = 0;
        }
      }
    } else {
      _simulation->_handlers->mpiFghExchangeStack[TDimension][TDirection] =
        FghMpiExchangeHandler::getSendHandler(
          _simulation->_memory.grid(),
          _simulation->_memory.parallelDistribution(),
          leftIndent,
          rightIndent);
      _simulation->_handlers->mpiPressureExchangeStack[TDimension][TDirection]
        = PressureMpiExchangeHandler::getReceiveHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.
        parallelDistribution(),
        leftIndent,
        rightIndent);

      for (int d = 0; d < Dimensions; ++d) {
        if (d != TDimension) {
          leftIndent(d) = 0;
        }
      }
    }

    _simulation->_handlers->mpiVelocityExchangeStack[TDimension][TDirection] =
      VelocityMpiExchangeHandler::getExchangeHandler(
        _simulation->_memory.grid(),
        _simulation->_memory.parallelDistribution(),
        leftIndent,
        rightIndent);
  }

private:
  Configuration* _configuration;
  SolverType*    _simulation;

  template <int TDimension>
  static TScalar*
  fghAccessor(CellAccessorType const& accessor) {
    return &accessor.fgh(TDimension);
  }

  static VectorDsType*
  velocityAccessor(CellAccessorType const& accessor) {
    return &accessor.velocity();
  }

  static ScalarType*
  pressureAccessor(CellAccessorType const& accessor) {
    return &accessor.pressure();
  }
};
}
}
