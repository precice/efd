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
template <typename Scalar, int TD>
class SimulationBuilder {
public:
  typedef FluidSimulation::FdSimulation<
      FluidSimulation::UniformGridGeometry<Scalar, TD>,
      StructuredMemory::IterableMemory<FluidSimulation::Cell<Scalar, TD>, TD>,
      Scalar,
      TD> Simulation;
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
      FluidSimulation::GhostLayer::Initialization::OutputFghAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      OutputFghInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, OutputFghInitializationAction, TD,
       TDimension, TDirection>
      OutputFghInitialization;
    typedef
      FluidSimulation::GhostLayer::PressureStencil::Handler
      <GridS, Scalar, TD, TDimension, TDirection>
      PressureStencil;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::MovingWallRhsAction
      <GridS, TD, TDimension, TDirection>
      MovingWallRhsAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, MovingWallRhsAction, TD,
       TDimension, TDirection>
      MovingWallRhsHandler;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::CopyPressureAction
      <GridS, TD, TDimension, TDirection>
      CopyPressureAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, CopyPressureAction, TD,
       TDimension, TDirection>
      CopyPressureHandler;
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
      FluidSimulation::GhostLayer::Initialization::OutputVelocityAction
      <typename GridS::Base, Scalar, TD, TDimension, TDirection>
      OutputVelocityInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, OutputVelocityInitializationAction, TD,
       TDimension, TDirection>
      OutputVelocityInitialization;
  };

public:
  SimulationBuilder(FluidSimulation::Configuration& parameters)
    : _parameters(parameters) {
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
    processorSize(0)  = parameters.parallel.numProcessors[0];
    globalCellSize(0) = parameters.geometry.sizeX;
    width(0)          = parameters.geometry.lengthX;
    processorSize(1)  = parameters.parallel.numProcessors[1];
    globalCellSize(1) = parameters.geometry.sizeY;
    width(1)          = parameters.geometry.lengthY;

    if (TD == 3) {
      processorSize(2)  = parameters.parallel.numProcessors[2];
      globalCellSize(2) = parameters.geometry.sizeZ;
      width(2)          = parameters.geometry.lengthZ;
    }

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    _parallelDistribution->initialize(rank,
                                      processorSize,
                                      globalCellSize);

    _simulation->_parameters.re()    = parameters.flow.Re;
    _simulation->_parameters.gamma() = parameters.solver.gamma;
    _simulation->_parameters.tau()   = parameters.timestep.tau;
    _simulation->_parameters.g(0)    = parameters.environment.gx;
    _simulation->_parameters.g(1)    = parameters.environment.gy;
    _simulation->_iterationLimit     = 200;
    _simulation->_timeLimit          = std::numeric_limits<Scalar>::max();

    if (TD == 3) {
      _simulation->_parameters.g(2) = parameters.environment.gz;
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
    typedef typename LeftHandlers::MovingWallFghInitializationAction
      FghAction;
    typedef typename LeftHandlers::MovingWallFghInitialization
      FghHandler;
    typedef typename LeftHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename LeftHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename LeftHandlers::MovingWallVelocityInitializationAction
      VelocityAction;
    typedef typename LeftHandlers::MovingWallVelocityInitialization
      VelocityHandler;
    typedef typename LeftHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[0][0] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[0][0],
        _parallelDistribution,
        new FghAction(_parameters));
    _handlers->rhsInitialization[0][0] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->velocityInitialization[0][0] =
      VelocityHandler::getHandler(
        &_grid->boundaries[0][0],
        _parallelDistribution,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->pressureStencilStack[0][0] =
      PressureStencil::getDirichletHandler(_grid, _parallelDistribution);
  }

  void
  _setLeftAsInput() {
    typedef GhostHandlers<0, 0> LeftHandlers;
    typedef typename LeftHandlers::InputFghInitializationAction
      FghAction;
    typedef typename LeftHandlers::InputFghInitialization
      FghHandler;
    typedef typename LeftHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename LeftHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename LeftHandlers::InputVelocityInitializationAction
      VelocityAction;
    typedef typename LeftHandlers::InputVelocityInitialization
      VelocityHandler;
    typedef typename LeftHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[0][0] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[0][0],
        _parallelDistribution,
        new FghAction(_parameters));
    _handlers->rhsInitialization[0][0] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->velocityInitialization[0][0] =
      VelocityHandler::getHandler(
        &_grid->boundaries[0][0],
        _parallelDistribution,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->pressureStencilStack[0][0] =
      PressureStencil::getDirichletHandler(_grid, _parallelDistribution);
  }

  void
  _setLeftAsOutput() {
    typedef GhostHandlers<0, 0> LeftHandlers;
    typedef typename LeftHandlers::OutputFghInitializationAction
      FghAction;
    typedef typename LeftHandlers::OutputFghInitialization
      FghHandler;
    typedef typename LeftHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename LeftHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename LeftHandlers::OutputVelocityInitializationAction
      VelocityAction;
    typedef typename LeftHandlers::OutputVelocityInitialization
      VelocityHandler;
    typedef typename LeftHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[0][0] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[0][0],
        _parallelDistribution,
        new FghAction());
    _handlers->rhsInitialization[0][0] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->velocityInitialization[0][0] =
      VelocityHandler::getHandler(
        &_grid->boundaries[0][0],
        _parallelDistribution,
        new VelocityAction(_maxVelocity));
    _handlers->pressureStencilStack[0][0] =
      PressureStencil::getNeumannHandler(_grid, _parallelDistribution);
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
    typedef typename RightHandlers::MovingWallFghInitializationAction
      FghAction;
    typedef typename RightHandlers::MovingWallFghInitialization
      FghHandler;
    typedef typename RightHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename RightHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename RightHandlers::CopyPressureAction
      CopyPressureAction;
    typedef typename RightHandlers::CopyPressureHandler
      CopyPressureHandler;
    typedef typename RightHandlers::MovingWallVelocityInitializationAction
      VelocityAction;
    typedef typename RightHandlers::MovingWallVelocityInitialization
      VelocityHandler;
    typedef typename RightHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[0][1] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[0][1],
        _parallelDistribution,
        new FghAction(_parameters));
    _handlers->rhsInitialization[0][1] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->pressureInitialization[0][1] =
      CopyPressureHandler::getHandler(
        _grid,
        new CopyPressureAction(_parallelDistribution));
    _handlers->velocityInitialization[0][1] =
      VelocityHandler::getHandler(
        &_grid->boundaries[0][1],
        _parallelDistribution,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->pressureStencilStack[0][1] =
      PressureStencil::getDirichletHandler(_grid, _parallelDistribution);
  }
  void
  _setRightAsInput() {
    typedef GhostHandlers<0, 1> RightHandlers;
    typedef typename RightHandlers::InputFghInitializationAction
      FghAction;
    typedef typename RightHandlers::InputFghInitialization
      FghHandler;
    typedef typename RightHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename RightHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename RightHandlers::CopyPressureAction
      CopyPressureAction;
    typedef typename RightHandlers::CopyPressureHandler
      CopyPressureHandler;
    typedef typename RightHandlers::InputVelocityInitializationAction
      VelocityAction;
    typedef typename RightHandlers::InputVelocityInitialization
      VelocityHandler;
    typedef typename RightHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[0][1] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[0][1],
        _parallelDistribution,
        new FghAction(_parameters));
    _handlers->rhsInitialization[0][1] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->pressureInitialization[0][1] =
      CopyPressureHandler::getHandler(
        _grid,
        new CopyPressureAction(_parallelDistribution));
    _handlers->velocityInitialization[0][1] =
      VelocityHandler::getHandler(
        &_grid->boundaries[0][1],
        _parallelDistribution,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->pressureStencilStack[0][1] =
      PressureStencil::getDirichletHandler(_grid, _parallelDistribution);
  }

  void
  _setRightAsOutput() {
    typedef GhostHandlers<0, 1> RightHandlers;
    typedef typename RightHandlers::OutputFghInitializationAction
      FghAction;
    typedef typename RightHandlers::OutputFghInitialization
      FghHandler;
    typedef typename RightHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename RightHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename RightHandlers::CopyPressureAction
      CopyPressureAction;
    typedef typename RightHandlers::CopyPressureHandler
      CopyPressureHandler;
    typedef typename RightHandlers::OutputVelocityInitializationAction
      VelocityAction;
    typedef typename RightHandlers::OutputVelocityInitialization
      VelocityHandler;
    typedef typename RightHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[0][1] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[0][1],
        _parallelDistribution,
        new FghAction());
    _handlers->rhsInitialization[0][1] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->pressureInitialization[0][1] =
      CopyPressureHandler::getHandler(
        _grid,
        new CopyPressureAction(_parallelDistribution));
    _handlers->velocityInitialization[0][1] =
      VelocityHandler::getHandler(
        &_grid->boundaries[0][1],
        _parallelDistribution,
        new VelocityAction(_maxVelocity));
    _handlers->pressureStencilStack[0][1] =
      PressureStencil::getNeumannHandler(_grid, _parallelDistribution);
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
    typedef typename BottomHandlers::MovingWallFghInitializationAction
      FghAction;
    typedef typename BottomHandlers::MovingWallFghInitialization
      FghHandler;
    typedef typename BottomHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename BottomHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename BottomHandlers::MovingWallVelocityInitializationAction
      VelocityAction;
    typedef typename BottomHandlers::MovingWallVelocityInitialization
      VelocityHandler;
    typedef typename BottomHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[1][0] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[1][0],
        _parallelDistribution,
        new FghAction(_parameters));
    _handlers->rhsInitialization[1][0] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->velocityInitialization[1][0] =
      VelocityHandler::getHandler(
        &_grid->boundaries[1][0],
        _parallelDistribution,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->pressureStencilStack[1][0] =
      PressureStencil::getDirichletHandler(_grid, _parallelDistribution);
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
    typedef typename TopHandlers::MovingWallFghInitializationAction
      FghAction;
    typedef typename TopHandlers::MovingWallFghInitialization
      FghHandler;
    typedef typename TopHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename TopHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename TopHandlers::CopyPressureAction
      CopyPressureAction;
    typedef typename TopHandlers::CopyPressureHandler
      CopyPressureHandler;
    typedef typename TopHandlers::MovingWallVelocityInitializationAction
      VelocityAction;
    typedef typename TopHandlers::MovingWallVelocityInitialization
      VelocityHandler;
    typedef typename TopHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[1][1] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[1][1],
        _parallelDistribution,
        new FghAction(_parameters));
    _handlers->rhsInitialization[1][1] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->pressureInitialization[1][1] =
      CopyPressureHandler::getHandler(
        _grid,
        new CopyPressureAction(_parallelDistribution));
    _handlers->velocityInitialization[1][1] =
      VelocityHandler::getHandler(
        &_grid->boundaries[1][1],
        _parallelDistribution,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->pressureStencilStack[1][1] =
      PressureStencil::getDirichletHandler(_grid, _parallelDistribution);
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
    typedef typename BackHandlers::MovingWallFghInitializationAction
      FghAction;
    typedef typename BackHandlers::MovingWallFghInitialization
      FghHandler;
    typedef typename BackHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename BackHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename BackHandlers::MovingWallVelocityInitializationAction
      VelocityAction;
    typedef typename BackHandlers::MovingWallVelocityInitialization
      VelocityHandler;
    typedef typename BackHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[2][0] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[2][0],
        _parallelDistribution,
        new FghAction(_parameters));
    _handlers->rhsInitialization[2][0] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->velocityInitialization[2][0] =
      VelocityHandler::getHandler(
        &_grid->boundaries[2][0],
        _parallelDistribution,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->pressureStencilStack[2][0] =
      PressureStencil::getDirichletHandler(_grid, _parallelDistribution);
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
    typedef typename FrontHandlers::MovingWallFghInitializationAction
      FghAction;
    typedef typename FrontHandlers::MovingWallFghInitialization
      FghHandler;
    typedef typename FrontHandlers::MovingWallRhsAction
      RhsAction;
    typedef typename FrontHandlers::MovingWallRhsHandler
      RhsHandler;
    typedef typename FrontHandlers::CopyPressureAction
      CopyPressureAction;
    typedef typename FrontHandlers::CopyPressureHandler
      CopyPressureHandler;
    typedef typename FrontHandlers::MovingWallVelocityInitializationAction
      VelocityAction;
    typedef typename FrontHandlers::MovingWallVelocityInitialization
      VelocityHandler;
    typedef typename FrontHandlers::PressureStencil
      PressureStencil;

    _handlers->fghInitialization[2][1] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[2][1],
        _parallelDistribution,
        new FghAction(_parameters));
    _handlers->rhsInitialization[2][1] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelDistribution));
    _handlers->pressureInitialization[2][1] =
      CopyPressureHandler::getHandler(
        _grid,
        new CopyPressureAction(_parallelDistribution));
    _handlers->velocityInitialization[2][1] =
      VelocityHandler::getHandler(
        &_grid->boundaries[2][1],
        _parallelDistribution,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->pressureStencilStack[2][1] =
      PressureStencil::getDirichletHandler(_grid, _parallelDistribution);
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

  FluidSimulation::Configuration& _parameters;
  Simulation*                     _simulation;
  GridS*                          _grid;
  ParallelTopologyS*              _parallelDistribution;
  GhostCellsHandlerS*             _handlers;
  VelocityS*                      _maxVelocity;
};
}
}
#endif
