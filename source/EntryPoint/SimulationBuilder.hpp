#ifndef FsiSimulation_EntryPoint_SimulationBuilder_hpp
#define FsiSimulation_EntryPoint_SimulationBuilder_hpp

#include "FluidSimulation/Cell.hpp"
#include "FluidSimulation/Configuration.hpp"
#include "FluidSimulation/FdSimulation.hpp"
#include "FluidSimulation/GhostLayer/InitializationActions.hpp"
#include "FluidSimulation/GhostLayer/PetscExchangeActions.hpp"
#include "FluidSimulation/GridGeometry.hpp"

#include "StructuredMemory/Accessor.hpp"
#include "StructuredMemory/Memory.hpp"

#include <Uni/StructuredGrid/Basic/MultiIndex>

namespace FsiSimulation {
namespace EntryPoint {
template <typename Scalar, int TD, int TSolverType = 0>
class SimulationBuilder {
public:
  typedef FluidSimulation::FdSimulation<
      FluidSimulation::UniformGridGeometry<Scalar, TD>,
      StructuredMemory::IterableMemory<FluidSimulation::Cell<Scalar, TD>, TD>,
      Scalar,
      TD,
      TSolverType> Simulation;
  typedef typename Simulation::GridType                 GridS;
  typedef typename Simulation::GridGeometryType         GridGeometryS;
  typedef typename Simulation::ParallelDistributionType ParallelTopologyS;
  typedef typename Simulation::GhostHandlersType        GhostCellsHandlerS;
  typedef typename Simulation::MemoryType               MemoryS;
  typedef typename Simulation::CellAccessorFactory
    CellAccessorFactoryS;
  typedef typename GridS::CellAccessorType CellAccessorS;
  typedef typename CellAccessorS::CellType CellS;
  typedef typename CellS::Velocity         VelocityS;
  typedef typename CellS::Pressure         PressureS;

  template <int TDimension, int TDirection>
  struct GhostHandlers {
    typedef
      FluidSimulation::GhostLayer::MpiExchange::Handler
      <Scalar, typename GridS::Base,
       SimulationBuilder::template fghAccessor<TDimension>,
       Scalar, TD, TDimension, TDirection>
      FghMpiExchange;
    typedef
      FluidSimulation::GhostLayer::MpiExchange::Handler
      <PressureS, typename GridS::Base,
       SimulationBuilder::pressureAccessor, Scalar, TD,
       TDimension, TDirection>
      PressureMpiExchange;
    typedef
      FluidSimulation::GhostLayer::MpiExchange::Handler
      <VelocityS, typename GridS::Base,
       SimulationBuilder::velocityAccessor, Scalar, TD,
       TDimension, TDirection>
      VelocityMpiExchange;

    typedef
      FluidSimulation::GhostLayer::Initialization::MovingWallFghAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      MovingWallFghInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, MovingWallFghInitializationAction, TD,
       TDimension, TDirection>
      MovingWallFghInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::InputFghAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      InputFghInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, InputFghInitializationAction, TD,
       TDimension, TDirection>
      InputFghInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::ParabolicInputFghAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      ParabolicInputFghInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, ParabolicInputFghInitializationAction, TD,
       TDimension, TDirection>
      ParabolicInputFghInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::OutputFghAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      OutputFghInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, OutputFghInitializationAction, TD,
       TDimension, TDirection>
      OutputFghInitialization;

    typedef
      FluidSimulation::GhostLayer::LsStencilGenerator::Handler
      <GridS, TDimension, TDirection>
      PpeStencilGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::ConstantRhsGenerationAction<
        CellAccessorS>
      PpeRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, PpeRhsGenerationAction, TDimension, TDirection>
      PpeRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::PpeRhsAcquiererAction
      PpeRhsAcquiererAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, PpeRhsAcquiererAction, TDimension, TDirection>
      PpeRhsAcquiererHandler;

    typedef
      FluidSimulation::GhostLayer::LsStencilGenerator::Handler
      <GridS, TDimension, TDirection>
      VpeStencilGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::ConstantRhsGenerationAction
      <CellAccessorS>
      VpeConstantRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VpeConstantRhsGenerationAction, TDimension, TDirection>
      VpeConstantRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeInputRhsGenerationAction
      <0, TDimension, TDirection>
      VxpeInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VxpeInputRhsGenerationAction, TDimension, TDirection>
      VxpeInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::
      VpeParabolicInputRhsGenerationAction
      <0, TDimension, TDirection>
      VxpeParabolicInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VxpeParabolicInputRhsGenerationAction, TDimension, TDirection>
      VxpeParabolicInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeRhsAcquiererAction<0>
      VxpeRhsAcquiererAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VxpeRhsAcquiererAction, TDimension, TDirection>
      VxpeRhsAcquiererHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeInputRhsGenerationAction
      <1, TDimension, TDirection>
      VypeInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VypeInputRhsGenerationAction, TDimension, TDirection>
      VypeInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::
      VpeParabolicInputRhsGenerationAction
      <1, TDimension, TDirection>
      VypeParabolicInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VypeParabolicInputRhsGenerationAction, TDimension, TDirection>
      VypeParabolicInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeRhsAcquiererAction<0>
      VypeRhsAcquiererAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VypeRhsAcquiererAction, TDimension, TDirection>
      VypeRhsAcquiererHandler;

    typedef
      FluidSimulation::GhostLayer::Initialization::MovingWallVelocityAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      MovingWallVelocityInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, MovingWallVelocityInitializationAction, TD,
       TDimension, TDirection>
      MovingWallVelocityInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::InputVelocityAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      InputVelocityInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, InputVelocityInitializationAction, TD,
       TDimension, TDirection>
      InputVelocityInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::ParabolicInputVelocityAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      ParabolicInputVelocityInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, ParabolicInputVelocityInitializationAction, TD,
       TDimension, TDirection>
      ParabolicInputVelocityInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::OutputVelocityAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      OutputVelocityInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, OutputVelocityInitializationAction, TD,
       TDimension, TDirection>
      OutputVelocityInitialization;

public:
    GhostHandlers(
      FluidSimulation::Configuration* configuration,
      GridS*                          grid,
      ParallelTopologyS*              parallelDistribution,
      GhostCellsHandlerS*             handlers,
      VelocityS*                      maxVelocity)
      : _configuration(configuration),
      _grid(grid),
      _parallelDistribution(parallelDistribution),
      _handlers(handlers),
      _maxVelocity(maxVelocity) {}

    inline void
    setAsInput() {
      _handlers->fghInitialization[TDimension][TDirection] =
        InputFghInitialization::getHandler(
          &_grid->indentedBoundaries[TDimension][TDirection],
          _parallelDistribution,
          new InputFghInitializationAction(_configuration));

      _handlers->ppeStencilGeneratorStack[TDimension][TDirection] =
        PpeStencilGenerationHandler::getNeumannMiddle(_grid,
                                                      _parallelDistribution);
      _handlers->ppeRhsGeneratorStack[TDimension][TDirection] =
        PpeRhsGenerationHandler::getHandler(_grid, _parallelDistribution,
                                            new PpeRhsGenerationAction(0.0));

      _handlers->vxpeRhsGeneratorStack[TDimension][TDirection] =
        VxpeInputRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VxpeInputRhsGenerationAction(_configuration));
      _handlers->vxpeRhsAcquiererStack[TDimension][TDirection] =
        VxpeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VxpeRhsAcquiererAction());

      _handlers->vypeRhsGeneratorStack[TDimension][TDirection] =
        VypeInputRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VypeInputRhsGenerationAction(_configuration));
      _handlers->vypeRhsAcquiererStack[TDimension][TDirection] =
        VypeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VypeRhsAcquiererAction());

      _handlers->velocityInitialization[TDimension][TDirection] =
        InputVelocityInitialization::getHandler(
          &_grid->boundaries[TDimension][TDirection],
          _parallelDistribution,
          new InputVelocityInitializationAction(_configuration, _maxVelocity));
    }

    inline void
    setAsParabolicInput() {
      _handlers->fghInitialization[TDimension][TDirection] =
        ParabolicInputFghInitialization::getHandler(
          &_grid->indentedBoundaries[TDimension][TDirection],
          _parallelDistribution,
          new ParabolicInputFghInitializationAction(_configuration));

      _handlers->ppeStencilGeneratorStack[TDimension][TDirection] =
        PpeStencilGenerationHandler::getNeumannMiddle(_grid,
                                                      _parallelDistribution);
      _handlers->ppeRhsGeneratorStack[TDimension][TDirection] =
        PpeRhsGenerationHandler::getHandler(_grid, _parallelDistribution,
                                            new PpeRhsGenerationAction(0.0));

      _handlers->vxpeRhsGeneratorStack[TDimension][TDirection] =
        VxpeParabolicInputRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VxpeParabolicInputRhsGenerationAction(_configuration));
      _handlers->vxpeRhsAcquiererStack[TDimension][TDirection] =
        VxpeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VxpeRhsAcquiererAction());

      _handlers->vypeRhsGeneratorStack[TDimension][TDirection] =
        VypeParabolicInputRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VypeParabolicInputRhsGenerationAction(_configuration));
      _handlers->vypeRhsAcquiererStack[TDimension][TDirection] =
        VypeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VypeRhsAcquiererAction());

      _handlers->velocityInitialization[TDimension][TDirection] =
        ParabolicInputVelocityInitialization::getHandler(
          &_grid->boundaries[TDimension][TDirection],
          _parallelDistribution,
          new ParabolicInputVelocityInitializationAction(_configuration,
                                                         _maxVelocity));
    }

    inline void
    setAsOutput() {
      _handlers->fghInitialization[TDimension][TDirection] =
        OutputFghInitialization::getHandler(
          &_grid->indentedBoundaries[TDimension][TDirection],
          _parallelDistribution,
          new OutputFghInitializationAction());

      _handlers->ppeStencilGeneratorStack[TDimension][TDirection] =
        PpeStencilGenerationHandler::getDirichletMiddle(_grid,
                                                        _parallelDistribution);
      _handlers->ppeRhsGeneratorStack[TDimension][TDirection] =
        PpeRhsGenerationHandler::getHandler(_grid, _parallelDistribution,
                                            new PpeRhsGenerationAction(0.0));

      _handlers->vxpeRhsGeneratorStack[TDimension][TDirection] =
        VpeConstantRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VpeConstantRhsGenerationAction(0.0));
      _handlers->vxpeRhsAcquiererStack[TDimension][TDirection] =
        VxpeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VxpeRhsAcquiererAction());

      _handlers->vypeRhsGeneratorStack[TDimension][TDirection] =
        VpeConstantRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VpeConstantRhsGenerationAction(0.0));
      _handlers->vypeRhsAcquiererStack[TDimension][TDirection] =
        VypeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VypeRhsAcquiererAction());

      _handlers->velocityInitialization[TDimension][TDirection] =
        OutputVelocityInitialization::getHandler(
          &_grid->boundaries[TDimension][TDirection],
          _parallelDistribution,
          new OutputVelocityInitializationAction(_maxVelocity));
    }

    inline void
    setVxpeInNonXDimensionAsInput() {
      _handlers->vxpeStencilGeneratorStack[TDimension][TDirection] =
        VpeStencilGenerationHandler::getDirichletMiddle(_grid,
                                                        _parallelDistribution);
    }

    inline void
    setVxpeInNonXDimensionAsParabolicInput() {
      _handlers->vxpeStencilGeneratorStack[TDimension][TDirection] =
        VpeStencilGenerationHandler::getDirichletMiddle(_grid,
                                                        _parallelDistribution);
    }

    inline void
    setVxpeInNonXDimensionAsOutput() {
      _handlers->vxpeStencilGeneratorStack[TDimension][TDirection] =
        VpeStencilGenerationHandler::getNeumannMiddle(_grid,
                                                      _parallelDistribution);
    }

    inline void
    setVypeInNonYDimensionAsInput() {
      _handlers->vypeStencilGeneratorStack[TDimension][TDirection] =
        VpeStencilGenerationHandler::getDirichletMiddle(_grid,
                                                        _parallelDistribution);
    }

    inline void
    setVypeInNonYDimensionAsParabolicInput() {
      _handlers->vypeStencilGeneratorStack[TDimension][TDirection] =
        VpeStencilGenerationHandler::getDirichletMiddle(_grid,
                                                        _parallelDistribution);
    }

    inline void
    setVypeInNonYDimensionAsOutput() {
      _handlers->vypeStencilGeneratorStack[TDimension][TDirection] =
        VpeStencilGenerationHandler::getNeumannMiddle(_grid,
                                                      _parallelDistribution);
    }

    inline void
    setPpeRhsAcquierer() {
      _handlers->ppeRhsAcquiererStack[TDimension][TDirection] =
        PpeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new PpeRhsAcquiererAction());
    }

private:
    FluidSimulation::Configuration* _configuration;
    GridS*                          _grid;
    ParallelTopologyS*              _parallelDistribution;
    GhostCellsHandlerS*             _handlers;
    VelocityS*                      _maxVelocity;
  };

public:
  SimulationBuilder(FluidSimulation::Configuration* configuration)
    : _configuration(configuration) {
    _simulation           = new Simulation();
    _grid                 = &_simulation->_grid;
    _parallelDistribution = &_simulation->_parallelDistribution;
    _handlers             = &_simulation->_ghostHandler;
    _maxVelocity          = &_simulation->_maxVelocity;

    typedef typename ParallelTopologyS::VectorDi VectorDi;
    typedef typename GridGeometryS::VectorDs     VectorDs;

    VectorDi processorSize;
    VectorDi globalCellSize;
    VectorDs width;
    processorSize(0)  = configuration->parallelizationSize(0);
    globalCellSize(0) = configuration->size(0);
    width(0)          = (Scalar)configuration->width(0);
    processorSize(1)  = configuration->parallelizationSize(1);
    globalCellSize(1) = configuration->size(1);
    width(1)          = (Scalar)configuration->width(1);

    if (TD == 3) {
      processorSize(2)  = configuration->parallelizationSize(2);
      globalCellSize(2) = configuration->size(2);
      width(2)          = (Scalar)configuration->width(2);
    }

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    _parallelDistribution->initialize(rank,
                                      processorSize,
                                      globalCellSize);

    _simulation->_parameters.re()    = configuration->re;
    _simulation->_parameters.gamma() = configuration->gamma;
    _simulation->_parameters.tau()   = configuration->tau;
    _simulation->_parameters.alpha() = configuration->alpha;
    _simulation->_parameters.immersedBoundaryMethod()
                                  = configuration->immersedBoundaryMethod;
    _simulation->_parameters.g(0) = configuration->environment(0);
    _simulation->_parameters.g(1) = configuration->environment(1);
    _simulation->_iterationLimit  = configuration->iterationLimit;
    _simulation->_timeLimit       = configuration->timeLimit;
    _simulation->_plotInterval    = configuration->plotInterval;

    if (TD == 3) {
      _simulation->_parameters.g(2) = configuration->environment(2);
    }

    VectorDi localSize(_parallelDistribution->localCellSize + 2 *
                       VectorDi::Ones());

    _simulation->setParameters(localSize, width);
  }

  Simulation*
  simulation() const {
    return _simulation;
  }

  void
  setLeftAsMoving() {
    if (_parallelDistribution->neighbors[0][0] < 0) {
      _setLeftAsMoving();
    } else {
      _setLeftAsMpiExchange();
    }
  }

  void
  setLeftAsInput() {
    if (_parallelDistribution->neighbors[0][0] < 0) {
      _setLeftAsInput();
    } else {
      _setLeftAsMpiExchange();
    }
  }

  void
  setLeftAsParabolicInput() {
    if (_parallelDistribution->neighbors[0][0] < 0) {
      _setLeftAsParabolicInput();
    } else {
      _setLeftAsMpiExchange();
    }
  }

  void
  setLeftAsOutput() {
    if (_parallelDistribution->neighbors[0][0] < 0) {
      _setLeftAsOutput();
    } else {
      _setLeftAsMpiExchange();
    }
  }

  void
  setRightAsMoving() {
    if (_parallelDistribution->neighbors[0][1] < 0) {
      _setRightAsMoving();
    } else {
      _setRightAsMpiExchange();
    }
  }

  void
  setRightAsInput() {
    if (_parallelDistribution->neighbors[0][1] < 0) {
      _setRightAsInput();
    } else {
      _setRightAsMpiExchange();
    }
  }

  void
  setRightAsOutput() {
    if (_parallelDistribution->neighbors[0][1] < 0) {
      _setRightAsOutput();
    } else {
      _setRightAsMpiExchange();
    }
  }

  void
  setBottomAsMoving() {
    if (_parallelDistribution->neighbors[1][0] < 0) {
      _setBottomAsMoving();
    } else {
      _setBottomAsMpiExchange();
    }
  }

  void
  setTopAsMoving() {
    if (_parallelDistribution->neighbors[1][1] < 0) {
      _setTopAsMoving();
    } else {
      _setTopAsMpiExchange();
    }
  }

  void
  setBackAsMoving() {
    if (_parallelDistribution->neighbors[2][0] < 0) {
      _setBackAsMoving();
    } else {
      _setBackAsMpiExchange();
    }
  }

  void
  setFrontAsMoving() {
    if (_parallelDistribution->neighbors[2][1] < 0) {
      _setFrontAsMoving();
    } else {
      setFrontAsMpiExchange();
    }
  }

private:
  void
  _setLeftAsMoving() {
    typedef GhostHandlers<0, 0> LeftHandlers;

    typedef typename LeftHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    LeftHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsInput();
    handlers.setVypeInNonYDimensionAsInput();

    _handlers->vxpeStencilGeneratorStack[0][0] =
      VpeStencilGenerationHandler::getDirichletLeft(_grid,
                                                    _parallelDistribution);
  }

  void
  _setLeftAsInput() {
    typedef GhostHandlers<0, 0> LeftHandlers;
    typedef typename LeftHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    LeftHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsInput();
    handlers.setVypeInNonYDimensionAsInput();

    _handlers->vxpeStencilGeneratorStack[0][0] =
      VpeStencilGenerationHandler::getDirichletLeft(_grid,
                                                    _parallelDistribution);
  }

  void
  _setLeftAsParabolicInput() {
    typedef GhostHandlers<0, 0> LeftHandlers;
    typedef typename LeftHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    LeftHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsParabolicInput();
    handlers.setVypeInNonYDimensionAsParabolicInput();

    _handlers->vxpeStencilGeneratorStack[0][0] =
      VpeStencilGenerationHandler::getDirichletLeft(_grid,
                                                    _parallelDistribution);
  }

  void
  _setLeftAsOutput() {
    typedef GhostHandlers<0, 0> LeftHandlers;

    typedef typename LeftHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    LeftHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsOutput();
    handlers.setVypeInNonYDimensionAsOutput();

    _handlers->vxpeStencilGeneratorStack[0][0] =
      VpeStencilGenerationHandler::getNeumannLeft(_grid,
                                                  _parallelDistribution);
  }

  void
  _setLeftAsMpiExchange() {
    typedef GhostHandlers<0, 0>                        LeftHandlers;
    typedef typename LeftHandlers::FghMpiExchange      FghHandler;
    typedef typename LeftHandlers::PressureMpiExchange PressureHandler;
    typedef typename LeftHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[0][0] =
      FghHandler::getReceiveHandler(_grid,
                                    _parallelDistribution,
                                    leftIndent,
                                    rightIndent);

    _handlers->mpiPressureExchangeStack[0][0] =
      PressureHandler::getSendHandler(_grid,
                                      _parallelDistribution,
                                      leftIndent,
                                      rightIndent);

    rightIndent(1) = 0;
    rightIndent(2) = 0;

    _handlers->mpiVelocityExchangeStack[0][0] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setRightAsMoving() {
    typedef GhostHandlers<0, 1> RightHandlers;

    typedef typename RightHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    RightHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsInput();
    handlers.setVypeInNonYDimensionAsInput();
    handlers.setPpeRhsAcquierer();

    _handlers->vxpeStencilGeneratorStack[0][1] =
      VpeStencilGenerationHandler::getDirichletRight(_grid,
                                                     _parallelDistribution);
  }

  void
  _setRightAsInput() {
    typedef GhostHandlers<0, 1> RightHandlers;

    typedef typename RightHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    RightHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsInput();
    handlers.setVypeInNonYDimensionAsInput();
    handlers.setPpeRhsAcquierer();

    _handlers->vxpeStencilGeneratorStack[0][1] =
      VpeStencilGenerationHandler::getDirichletRight(_grid,
                                                     _parallelDistribution);
  }

  void
  _setRightAsOutput() {
    typedef GhostHandlers<0, 1> RightHandlers;

    typedef typename RightHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    RightHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsOutput();
    handlers.setVypeInNonYDimensionAsOutput();
    handlers.setPpeRhsAcquierer();

    _handlers->vxpeStencilGeneratorStack[0][1] =
      VpeStencilGenerationHandler::getNeumannRight(_grid,
                                                   _parallelDistribution);
  }

  void
  _setRightAsMpiExchange() {
    typedef GhostHandlers<0, 1>                         RightHandlers;
    typedef typename RightHandlers::FghMpiExchange      FghHandler;
    typedef typename RightHandlers::PressureMpiExchange PressureHandler;
    typedef typename RightHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[0][1] =
      FghHandler::getSendHandler(_grid,
                                 _parallelDistribution,
                                 leftIndent,
                                 rightIndent);

    _handlers->mpiPressureExchangeStack[0][1] =
      PressureHandler::getReceiveHandler(_grid,
                                         _parallelDistribution,
                                         leftIndent,
                                         rightIndent);

    leftIndent(1) = 0;
    leftIndent(2) = 0;

    _handlers->mpiVelocityExchangeStack[0][1] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setBottomAsMoving() {
    typedef GhostHandlers<1, 0> BottomHandlers;
    typedef typename BottomHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    BottomHandlers handlers(_configuration,
                            _grid,
                            _parallelDistribution,
                            _handlers,
                            _maxVelocity);

    handlers.setAsInput();
    handlers.setVxpeInNonXDimensionAsInput();

    _handlers->vypeStencilGeneratorStack[1][0] =
      VpeStencilGenerationHandler::getDirichletLeft(_grid,
                                                    _parallelDistribution);
  }

  void
  _setBottomAsMpiExchange() {
    typedef GhostHandlers<1, 0>                          BottomHandlers;
    typedef typename BottomHandlers::FghMpiExchange      FghHandler;
    typedef typename BottomHandlers::PressureMpiExchange PressureHandler;
    typedef typename BottomHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[1][0] =
      FghHandler::getReceiveHandler(_grid,
                                    _parallelDistribution,
                                    leftIndent,
                                    rightIndent);

    _handlers->mpiPressureExchangeStack[1][0] =
      PressureHandler::getSendHandler(_grid,
                                      _parallelDistribution,
                                      leftIndent,
                                      rightIndent);

    rightIndent(0) = 0;
    rightIndent(2) = 0;

    _handlers->mpiVelocityExchangeStack[1][0] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setTopAsMoving() {
    typedef GhostHandlers<1, 1> TopHandlers;

    typedef typename TopHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    TopHandlers handlers(_configuration,
                         _grid,
                         _parallelDistribution,
                         _handlers,
                         _maxVelocity);

    handlers.setAsInput();
    handlers.setVxpeInNonXDimensionAsInput();
    handlers.setPpeRhsAcquierer();

    _handlers->vypeStencilGeneratorStack[1][1] =
      VpeStencilGenerationHandler::getDirichletRight(_grid,
                                                     _parallelDistribution);
  }

  void
  _setTopAsMpiExchange() {
    typedef GhostHandlers<1, 1>                       TopHandlers;
    typedef typename TopHandlers::FghMpiExchange      FghHandler;
    typedef typename TopHandlers::PressureMpiExchange PressureHandler;
    typedef typename TopHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[1][1] =
      FghHandler::getSendHandler(_grid,
                                 _parallelDistribution,
                                 leftIndent,
                                 rightIndent);

    _handlers->mpiPressureExchangeStack[1][1] =
      PressureHandler::getReceiveHandler(_grid,
                                         _parallelDistribution,
                                         leftIndent,
                                         rightIndent);

    leftIndent(0) = 0;
    leftIndent(2) = 0;

    _handlers->mpiVelocityExchangeStack[1][1] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setBackAsMoving() {
    typedef GhostHandlers<2, 0> BackHandlers;
    typedef typename BackHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    BackHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsInput();
    handlers.setVxpeInNonXDimensionAsInput();
    handlers.setVypeInNonYDimensionAsInput();

    // _handlers->vxpeStencilGeneratorStack[2][0] =
    //   VpeStencilGenerationHandler::getDirichletLeft(_grid,
    //                                                 _parallelDistribution);
  }

  void
  _setBackAsMpiExchange() {
    typedef GhostHandlers<2, 0>                        BackHandlers;
    typedef typename BackHandlers::FghMpiExchange      FghHandler;
    typedef typename BackHandlers::PressureMpiExchange PressureHandler;
    typedef typename BackHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[2][0] =
      FghHandler::getReceiveHandler(_grid,
                                    _parallelDistribution,
                                    leftIndent,
                                    rightIndent);

    _handlers->mpiPressureExchangeStack[2][0] =
      PressureHandler::getSendHandler(_grid,
                                      _parallelDistribution,
                                      leftIndent,
                                      rightIndent);

    rightIndent(0) = 0;
    rightIndent(1) = 0;

    _handlers->mpiVelocityExchangeStack[2][0] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setFrontAsMoving() {
    typedef GhostHandlers<2, 1> FrontHandlers;

    typedef typename FrontHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    FrontHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsInput();
    handlers.setVxpeInNonXDimensionAsInput();
    handlers.setVypeInNonYDimensionAsInput();
    handlers.setPpeRhsAcquierer();

    // _handlers->vxpeStencilGeneratorStack[2][1] =
    //   VpeStencilGenerationHandler::getDirichletRight(_grid,
    //                                                  _parallelDistribution);
  }

  void
  setFrontAsMpiExchange() {
    typedef GhostHandlers<2, 1>                         FrontHandlers;
    typedef typename FrontHandlers::FghMpiExchange      FghHandler;
    typedef typename FrontHandlers::PressureMpiExchange PressureHandler;
    typedef typename FrontHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[2][1] =
      FghHandler::getSendHandler(_grid,
                                 _parallelDistribution,
                                 leftIndent,
                                 rightIndent);

    _handlers->mpiPressureExchangeStack[2][1] =
      PressureHandler::getReceiveHandler(_grid,
                                         _parallelDistribution,
                                         leftIndent,
                                         rightIndent);

    leftIndent(0) = 0;
    leftIndent(1) = 0;

    _handlers->mpiVelocityExchangeStack[2][1] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  template <int TDimension>
  static Scalar*
  fghAccessor(CellAccessorS const& accessor) {
    return &accessor.currentCell()->fgh(TDimension);
  }

  static VelocityS*
  velocityAccessor(CellAccessorS const& accessor) {
    return &accessor.currentCell()->velocity();
  }

  static PressureS*
  pressureAccessor(CellAccessorS const& accessor) {
    return &accessor.currentCell()->pressure();
  }

  FluidSimulation::Configuration* _configuration;
  Simulation*                     _simulation;
  GridS*                          _grid;
  ParallelTopologyS*              _parallelDistribution;
  GhostCellsHandlerS*             _handlers;
  VelocityS*                      _maxVelocity;
};
}
}
#endif
