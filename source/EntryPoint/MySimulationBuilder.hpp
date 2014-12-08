#ifndef FsiSimulation_EntryPoint_MySimulationBuilder_hpp
#define FsiSimulation_EntryPoint_MySimulationBuilder_hpp

#include "Cell.hpp"
#include "GridGeometry.hpp"
#include "MyTemplateSimulation.hpp"
#include "Solvers/GhostInitializationActions.hpp"
#include "Solvers/GhostPetscExchangeActions.hpp"
#include "StructuredMemory/Accessor.hpp"
#include "StructuredMemory/Memory.hpp"

#include <Uni/StructuredGrid/Basic/MultiIndex>

namespace FsiSimulation {
namespace EntryPoint {
template <typename Scalar, int D>
class MySimulationBuilder {
public:
  typedef MyTemplateSimulation<
      UniformGridGeometry<Scalar, D>,
      StructuredMemory::IterableMemory<Cell<Scalar, D>, D>,
      Scalar,
      D> Simulation;
  typedef typename Simulation::SpecializedGrid              GridS;
  typedef typename Simulation::GridGeometry                 GridGeometryS;
  typedef typename Simulation::SpecializedParallelTopology  ParallelTopologyS;
  typedef typename Simulation::SpecializedGhostCellsHandler GhostCellsHandlerS;
  typedef typename Simulation::Memory                       MemoryS;
  typedef typename Simulation::CellAccessorFactory
    CellAccessorFactoryS;
  typedef typename GridS::CellAccessor CellAccessorS;
  typedef typename CellAccessorS::Cell CellS;
  typedef typename CellS::Velocity     VelocityS;
  typedef typename CellS::Pressure     PressureS;

  template <int TDimension, int TDirection>
  struct GhostHandlers {
    typedef
      Solvers::Ghost::MpiExchange::Handler
      <Scalar, typename GridS::Base,
       MySimulationBuilder::template fghAccessor<TDimension>,
       Scalar, D, TDimension, TDirection>
      FghMpiExchange;
    typedef
      Solvers::Ghost::MpiExchange::Handler
      <PressureS, typename GridS::Base,
       MySimulationBuilder::pressureAccessor, Scalar, D,
       TDimension, TDirection>
      PressureMpiExchange;
    typedef
      Solvers::Ghost::MpiExchange::Handler
      <VelocityS, typename GridS::Base,
       MySimulationBuilder::velocityAccessor, Scalar, D,
       TDimension, TDirection>
      VelocityMpiExchange;
    typedef
      Solvers::Ghost::Initialization::MovingWallFghAction
      <typename GridS::Base, Scalar, D, TDimension, TDirection>
      MovingWallFghInitializationAction;
    typedef
      Solvers::Ghost::Initialization::Handler
      <typename GridS::Base, MovingWallFghInitializationAction, D,
       TDimension, TDirection>
      MovingWallFghInitialization;
    typedef
      Solvers::Ghost::PetscExchange::MovingWallRhsAction
      <GridS, D, TDimension, TDirection>
      MovingWallRhsAction;
    typedef
      Solvers::Ghost::PetscExchange::Handler
      <GridS, MovingWallRhsAction, D,
       TDimension, TDirection>
      MovingWallRhsHandler;
    typedef
      Solvers::Ghost::PetscExchange::CopyPressureAction
      <GridS, D, TDimension, TDirection>
      CopyPressureAction;
    typedef
      Solvers::Ghost::PetscExchange::Handler
      <GridS, CopyPressureAction, D,
       TDimension, TDirection>
      CopyPressureHandler;
    typedef
      Solvers::Ghost::Initialization::MovingWallVelocityAction
      <typename GridS::Base, Scalar, D, TDimension, TDirection>
      MovingWallVelocityInitializationAction;
    typedef
      Solvers::Ghost::Initialization::Handler
      <typename GridS::Base, MovingWallVelocityInitializationAction, D,
       TDimension, TDirection>
      MovingWallVelocityInitialization;
    typedef
      Solvers::Ghost::PressureStencil::Handler
      <GridS, Scalar, D, TDimension, TDirection>
      PressureStencil;
  };

public:
  MySimulationBuilder(Parameters& parameters)
    : _parameters(parameters) {
    _simulation       = new Simulation();
    _grid             = &_simulation->_grid;
    _parallelTopology = &_simulation->_parallelTopology;
    _handlers         = &_simulation->_ghostCellsHandler;
    _maxVelocity      = &_simulation->_maxVelocity;

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

    if (D == 3) {
      processorSize(2)  = parameters.parallel.numProcessors[2];
      globalCellSize(2) = parameters.geometry.sizeZ;
      width(2)          = parameters.geometry.lengthZ;
    }

    _parallelTopology->initialize(parameters.parallel.rank,
                                  processorSize,
                                  globalCellSize);

    _simulation->_parameters.re()    = parameters.flow.Re;
    _simulation->_parameters.gamma() = parameters.solver.gamma;
    _simulation->_parameters.tau()   = parameters.timestep.tau;
    _simulation->_parameters.g(0)    = parameters.environment.gx;
    _simulation->_parameters.g(1)    = parameters.environment.gy;
    _simulation->_iterationLimit     = 200;
    _simulation->_timeLimit          = std::numeric_limits<Scalar>::max();

    if (D == 3) {
      _simulation->_parameters.g(2) = parameters.environment.gz;
    }

    VectorDi localSize(_parallelTopology->localSize + 2 * VectorDi::Ones());

    _simulation->setParameters(localSize, width);
  }

  Simulation*
  simulation() const {
    return _simulation;
  }

  void
  setLeftWallMoving() {
    if (_parallelTopology->neighbors[0][0] < 0) {
      _setLeftWallMoving();
    } else {
      _setLeftMpiExchange();
    }
  }

  void
  setRightWallMoving() {
    if (_parallelTopology->neighbors[0][1] < 0) {
      _setRightWallMoving();
    } else {
      _setRightMpiExchange();
    }
  }

  void
  setBottomWallMoving() {
    if (_parallelTopology->neighbors[1][0] < 0) {
      _setBottomWallMoving();
    } else {
      _setBottomMpiExchange();
    }
  }

  void
  setTopWallMoving() {
    if (_parallelTopology->neighbors[1][1] < 0) {
      _setTopWallMoving();
    } else {
      _setTopMpiExchange();
    }
  }

  void
  setBackWallMoving() {
    if (_parallelTopology->neighbors[2][0] < 0) {
      _setBackWallMoving();
    } else {
      _setBackMpiExchange();
    }
  }

  void
  setFrontWallMoving() {
    if (_parallelTopology->neighbors[2][1] < 0) {
      _setFrontWallMoving();
    } else {
      _setFrontMpiExchange();
    }
  }

private:
  void
  _setLeftWallMoving() {
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

    _handlers->_fghInitialization[0][0] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[0][0],
        _parallelTopology,
        new FghAction(_parameters));
    _handlers->_rhsInitialization[0][0] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelTopology));
    _handlers->_velocityInitialization[0][0] =
      VelocityHandler::getHandler(
        &_grid->boundaries[0][0],
        _parallelTopology,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->_pressureStencilStack[0][0] =
      PressureStencil::getDirichletHandler(_grid, _parallelTopology);
  }

  void
  _setLeftMpiExchange() {
    typedef GhostHandlers<0, 0>                        LeftHandlers;
    typedef typename LeftHandlers::FghMpiExchange      FghHandler;
    typedef typename LeftHandlers::PressureMpiExchange PressureHandler;
    typedef typename LeftHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->_mpiFghExchangeStack[0][0] =
      FghHandler::getReceiveHandler(_grid,
                                    _parallelTopology,
                                    leftIndent,
                                    rightIndent);

    _handlers->_mpiPressureExchangeStack[0][0] =
      PressureHandler::getSendHandler(_grid,
                                      _parallelTopology,
                                      leftIndent,
                                      rightIndent);

    rightIndent(1) = 0;
    rightIndent(2) = 0;

    _handlers->_mpiVelocityExchangeStack[0][0] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelTopology,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setRightWallMoving() {
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

    _handlers->_fghInitialization[0][1] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[0][1],
        _parallelTopology,
        new FghAction(_parameters));
    _handlers->_rhsInitialization[0][1] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelTopology));
    _handlers->_pressureInitialization[0][1] =
      CopyPressureHandler::getHandler(
        _grid,
        new CopyPressureAction(_parallelTopology));
    _handlers->_velocityInitialization[0][1] =
      VelocityHandler::getHandler(
        &_grid->boundaries[0][1],
        _parallelTopology,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->_pressureStencilStack[0][1] =
      PressureStencil::getDirichletHandler(_grid, _parallelTopology);
  }

  void
  _setRightMpiExchange() {
    typedef GhostHandlers<0, 1>                         RightHandlers;
    typedef typename RightHandlers::FghMpiExchange      FghHandler;
    typedef typename RightHandlers::PressureMpiExchange PressureHandler;
    typedef typename RightHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->_mpiFghExchangeStack[0][1] =
      FghHandler::getSendHandler(_grid,
                                 _parallelTopology,
                                 leftIndent,
                                 rightIndent);

    _handlers->_mpiPressureExchangeStack[0][1] =
      PressureHandler::getReceiveHandler(_grid,
                                         _parallelTopology,
                                         leftIndent,
                                         rightIndent);

    leftIndent(1) = 0;
    leftIndent(2) = 0;

    _handlers->_mpiVelocityExchangeStack[0][1] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelTopology,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setBottomWallMoving() {
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

    _handlers->_fghInitialization[1][0] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[1][0],
        _parallelTopology,
        new FghAction(_parameters));
    _handlers->_rhsInitialization[1][0] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelTopology));
    _handlers->_velocityInitialization[1][0] =
      VelocityHandler::getHandler(
        &_grid->boundaries[1][0],
        _parallelTopology,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->_pressureStencilStack[1][0] =
      PressureStencil::getDirichletHandler(_grid, _parallelTopology);
  }

  void
  _setBottomMpiExchange() {
    typedef GhostHandlers<1, 0>                          BottomHandlers;
    typedef typename BottomHandlers::FghMpiExchange      FghHandler;
    typedef typename BottomHandlers::PressureMpiExchange PressureHandler;
    typedef typename BottomHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->_mpiFghExchangeStack[1][0] =
      FghHandler::getReceiveHandler(_grid,
                                    _parallelTopology,
                                    leftIndent,
                                    rightIndent);

    _handlers->_mpiPressureExchangeStack[1][0] =
      PressureHandler::getSendHandler(_grid,
                                      _parallelTopology,
                                      leftIndent,
                                      rightIndent);

    rightIndent(0) = 0;
    rightIndent(2) = 0;

    _handlers->_mpiVelocityExchangeStack[1][0] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelTopology,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setTopWallMoving() {
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

    _handlers->_fghInitialization[1][1] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[1][1],
        _parallelTopology,
        new FghAction(_parameters));
    _handlers->_rhsInitialization[1][1] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelTopology));
    _handlers->_pressureInitialization[1][1] =
      CopyPressureHandler::getHandler(
        _grid,
        new CopyPressureAction(_parallelTopology));
    _handlers->_velocityInitialization[1][1] =
      VelocityHandler::getHandler(
        &_grid->boundaries[1][1],
        _parallelTopology,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->_pressureStencilStack[1][1] =
      PressureStencil::getDirichletHandler(_grid, _parallelTopology);
  }

  void
  _setTopMpiExchange() {
    typedef GhostHandlers<1, 1>                       TopHandlers;
    typedef typename TopHandlers::FghMpiExchange      FghHandler;
    typedef typename TopHandlers::PressureMpiExchange PressureHandler;
    typedef typename TopHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->_mpiFghExchangeStack[1][1] =
      FghHandler::getSendHandler(_grid,
                                 _parallelTopology,
                                 leftIndent,
                                 rightIndent);

    _handlers->_mpiPressureExchangeStack[1][1] =
      PressureHandler::getReceiveHandler(_grid,
                                         _parallelTopology,
                                         leftIndent,
                                         rightIndent);

    leftIndent(0) = 0;
    leftIndent(2) = 0;

    _handlers->_mpiVelocityExchangeStack[1][1] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelTopology,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setBackWallMoving() {
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

    _handlers->_fghInitialization[2][0] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[2][0],
        _parallelTopology,
        new FghAction(_parameters));
    _handlers->_rhsInitialization[2][0] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelTopology));
    _handlers->_velocityInitialization[2][0] =
      VelocityHandler::getHandler(
        &_grid->boundaries[2][0],
        _parallelTopology,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->_pressureStencilStack[2][0] =
      PressureStencil::getDirichletHandler(_grid, _parallelTopology);
  }

  void
  _setBackMpiExchange() {
    typedef GhostHandlers<2, 0>                        BackHandlers;
    typedef typename BackHandlers::FghMpiExchange      FghHandler;
    typedef typename BackHandlers::PressureMpiExchange PressureHandler;
    typedef typename BackHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->_mpiFghExchangeStack[2][0] =
      FghHandler::getReceiveHandler(_grid,
                                    _parallelTopology,
                                    leftIndent,
                                    rightIndent);

    _handlers->_mpiPressureExchangeStack[2][0] =
      PressureHandler::getSendHandler(_grid,
                                      _parallelTopology,
                                      leftIndent,
                                      rightIndent);

    rightIndent(0) = 0;
    rightIndent(1) = 0;

    _handlers->_mpiVelocityExchangeStack[2][0] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelTopology,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setFrontWallMoving() {
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

    _handlers->_fghInitialization[2][1] =
      FghHandler::getHandler(
        &_grid->indentedBoundaries[2][1],
        _parallelTopology,
        new FghAction(_parameters));
    _handlers->_rhsInitialization[2][1] =
      RhsHandler::getHandler(
        _grid,
        new RhsAction(_parameters, _parallelTopology));
    _handlers->_pressureInitialization[2][1] =
      CopyPressureHandler::getHandler(
        _grid,
        new CopyPressureAction(_parallelTopology));
    _handlers->_velocityInitialization[2][1] =
      VelocityHandler::getHandler(
        &_grid->boundaries[2][1],
        _parallelTopology,
        new VelocityAction(_parameters, _maxVelocity));
    _handlers->_pressureStencilStack[2][1] =
      PressureStencil::getDirichletHandler(_grid, _parallelTopology);
  }

  void
  _setFrontMpiExchange() {
    typedef GhostHandlers<2, 1>                         FrontHandlers;
    typedef typename FrontHandlers::FghMpiExchange      FghHandler;
    typedef typename FrontHandlers::PressureMpiExchange PressureHandler;
    typedef typename FrontHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->_mpiFghExchangeStack[2][1] =
      FghHandler::getSendHandler(_grid,
                                 _parallelTopology,
                                 leftIndent,
                                 rightIndent);

    _handlers->_mpiPressureExchangeStack[2][1] =
      PressureHandler::getReceiveHandler(_grid,
                                         _parallelTopology,
                                         leftIndent,
                                         rightIndent);

    leftIndent(0) = 0;
    leftIndent(1) = 0;

    _handlers->_mpiVelocityExchangeStack[2][1] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelTopology,
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

  Parameters&         _parameters;
  Simulation*         _simulation;
  GridS*              _grid;
  ParallelTopologyS*  _parallelTopology;
  GhostCellsHandlerS* _handlers;
  VelocityS*          _maxVelocity;
};
}
}
#endif
