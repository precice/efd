#include "SimulationFactory.hpp"

#include "Cell.hpp"
#include "GridGeometry.hpp"
#include "MyTemplateSimulation.hpp"
#include "Solvers/GhostInitializationActions.hpp"
#include "StructuredMemory/Accessor.hpp"
#include "StructuredMemory/Memory.hpp"

#include <Uni/StructuredGrid/Basic/MultiIndex>
//

using FsiSimulation::EntryPoint::SimulationFactory;

namespace FsiSimulation {
namespace EntryPoint {
namespace Private {
template <typename Scalar, int D>
using
  SpecializedSimulationT =
    MyTemplateSimulation<
      UniformGridGeometry<Scalar, D>,
      StructuredMemory::IterableMemory<Cell<Scalar, D>, D>,
      Scalar,
      D>;
template <typename Scalar, int D>
using GridT = typename SpecializedSimulationT<Scalar, D>::SpecializedGrid;
template <typename Scalar, int D>
using VectorDiT =
        typename SpecializedSimulationT<Scalar, D>::
        SpecializedParallelTopology::VectorDi;
template <typename Scalar, int D>
using VectorDsT =
        typename SpecializedSimulationT<Scalar, D>::GridGeometry::
        VectorDs;

template <typename Scalar, int D>
typename GridT<Scalar, D>::CellAccessor::Cell::Velocity *
velocityAccessor(typename GridT<Scalar, D>::CellAccessor const & accessor) {
  return &accessor.currentCell()->velocity();
}

template <typename Scalar, int D>
typename GridT<Scalar, D>::CellAccessor::Cell::Pressure *
pressureAccessor(typename GridT<Scalar, D>::CellAccessor const & accessor) {
  return &accessor.currentCell()->pressure();
}

template <typename Scalar, int D>
SimulationFactory::Simulation*
createUniformGridFromTemplate(Parameters& parameters) {
  typedef SpecializedSimulationT<Scalar, D>
    SpecializedSimulation;
  typedef typename SpecializedSimulation::SpecializedGrid
    Grid;
  typedef typename SpecializedSimulation::SpecializedParallelTopology::VectorDi
    VectorDi;
  typedef typename SpecializedSimulation::GridGeometry::VectorDs
    VectorDs;

  auto simulation = new SpecializedSimulation();

  parameters.walls.velocities[0][0](0) = parameters.walls.vectorLeft[0];
  parameters.walls.velocities[0][0](1) = parameters.walls.vectorLeft[1];
  parameters.walls.velocities[0][0](2) = parameters.walls.vectorLeft[2];

  parameters.walls.velocities[0][1](0) = parameters.walls.vectorRight[0];
  parameters.walls.velocities[0][1](1) = parameters.walls.vectorRight[1];
  parameters.walls.velocities[0][1](2) = parameters.walls.vectorRight[2];

  parameters.walls.velocities[1][0](0) = parameters.walls.vectorBottom[0];
  parameters.walls.velocities[1][0](1) = parameters.walls.vectorBottom[1];
  parameters.walls.velocities[1][0](2) = parameters.walls.vectorBottom[2];

  parameters.walls.velocities[1][1](0) = parameters.walls.vectorTop[0];
  parameters.walls.velocities[1][1](1) = parameters.walls.vectorTop[1];
  parameters.walls.velocities[1][1](2) = parameters.walls.vectorTop[2];

  parameters.walls.velocities[2][0](0) = parameters.walls.vectorBack[0];
  parameters.walls.velocities[2][0](1) = parameters.walls.vectorBack[1];
  parameters.walls.velocities[2][0](2) = parameters.walls.vectorBack[2];

  parameters.walls.velocities[2][1](0) = parameters.walls.vectorFront[0];
  parameters.walls.velocities[2][1](1) = parameters.walls.vectorFront[1];
  parameters.walls.velocities[2][1](2) = parameters.walls.vectorFront[2];

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

  simulation->_parallelTopology.initialize(parameters.parallel.rank,
                                           processorSize,
                                           globalCellSize);

  simulation->_parameters.re()    = parameters.flow.Re;
  simulation->_parameters.gamma() = parameters.solver.gamma;
  simulation->_parameters.tau()   = parameters.timestep.tau;
  simulation->_parameters.g(0)    = parameters.environment.gx;
  simulation->_parameters.g(1)    = parameters.environment.gy;
  simulation->_iterationLimit     = 100;
  simulation->_timeLimit          = std::numeric_limits<Scalar>::max();

  if (D == 3) {
    simulation->_parameters.g(2) = parameters.environment.gz;
  }

  typedef
    Solvers::Ghost::PressureStencil::DirichletStack
    <Grid, Scalar, D> DirichletStack;
  auto dirichletStack = DirichletStack::create(
    &simulation->_grid,
    &simulation->_parallelTopology);

  for (int d = 0; d < D; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      simulation->_ghostCellsHandler._pressureStencilStack[d][d2] =
        dirichletStack[d][d2];
    }
  }

  typedef
    Solvers::Ghost::MpiExchange::Stack
    < typename Grid::CellAccessor::Cell::Velocity,
    Grid,
    Scalar,
    D,
    velocityAccessor < Scalar, D >>
    MpiVelocityExchangeStack;

  typedef
    Solvers::Ghost::MpiExchange::Stack
    < typename Grid::CellAccessor::Cell::Pressure,
    Grid,
    Scalar,
    D,
    pressureAccessor < Scalar, D >>
    MpiPressureExchangeStack;

  auto mpiVelocityExchangeStack =
    MpiVelocityExchangeStack::create(
      &simulation->_grid,
      &simulation->_parallelTopology);

  auto mpiPressureExchangeStack =
    MpiPressureExchangeStack::create(
      &simulation->_grid,
      &simulation->_parallelTopology);

  for (int d = 0; d < D; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      if (simulation->_parallelTopology.neighbors[d][d2] >= 0) {
        simulation->_ghostCellsHandler._mpiVelocityExchangeStack[d][d2] =
          mpiVelocityExchangeStack[d][d2];
        simulation->_ghostCellsHandler._mpiPressureExchangeStack[d][d2] =
          mpiPressureExchangeStack[d][d2];
      }
    }
  }

  simulation->_gridGeometry.initialize(width,
                                       simulation->_parallelTopology.
                                       globalSize,
                                       simulation->_parallelTopology.corner);

  typedef
    Solvers::Ghost::Initialization::MovingWallFghAction
    <Grid, Scalar, D, 0, 0>
    LeftMovingWallFghAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, LeftMovingWallFghAction, D, 0, 0>
    InitializeLeftFghHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallFghAction
    <Grid, Scalar, D, 0, 1>
    RightMovingWallFghAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, RightMovingWallFghAction, D, 0, 1>
    InitializeRightFghHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallFghAction
    <Grid, Scalar, D, 1, 0>
    BottomMovingWallFghAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, BottomMovingWallFghAction, D, 1, 0>
    InitializeBottomFghHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallFghAction
    <Grid, Scalar, D, 1, 1>
    TopMovingWallFghAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, TopMovingWallFghAction, D, 1, 1>
    InitializeTopFghHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallFghAction
    <Grid, Scalar, D, 2, 0>
    BackMovingWallFghAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, BackMovingWallFghAction, D, 2, 0>
    InitializeBackFghHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallFghAction
    <Grid, Scalar, D, 2, 1>
    FrontMovingWallFghAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, FrontMovingWallFghAction, D, 2, 1>
    InitializeFrontFghHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallVelocityAction
    <Grid, Scalar, D, 0, 0>
    LeftMovingWallVelocityAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, LeftMovingWallVelocityAction, D, 0, 0>
    InitializeLeftVelocityHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallVelocityAction
    <Grid, Scalar, D, 0, 1>
    RightMovingWallVelocityAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, RightMovingWallVelocityAction, D, 0, 1>
    InitializeRightVelocityHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallVelocityAction
    <Grid, Scalar, D, 1, 0>
    BottomMovingWallVelocityAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, BottomMovingWallVelocityAction, D, 1, 0>
    InitializeBottomVelocityHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallVelocityAction
    <Grid, Scalar, D, 1, 1>
    TopMovingWallVelocityAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, TopMovingWallVelocityAction, D, 1, 1>
    InitializeTopVelocityHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallVelocityAction
    <Grid, Scalar, D, 2, 0>
    BackMovingWallVelocityAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, BackMovingWallVelocityAction, D, 2, 0>
    InitializeBackVelocityHandler;

  typedef
    Solvers::Ghost::Initialization::MovingWallVelocityAction
    <Grid, Scalar, D, 2, 1>
    FrontMovingWallVelocityAction;
  typedef
    Solvers::Ghost::Initialization::Handler
    <Grid, Scalar, FrontMovingWallVelocityAction, D, 2, 1>
    InitializeFrontVelocityHandler;

  simulation->_ghostCellsHandler._fghInitialization[0][0]
    = InitializeLeftFghHandler::getHandler(
    &simulation->_grid,
    &simulation->_parallelTopology,
    new LeftMovingWallFghAction(parameters));

  simulation->_ghostCellsHandler._fghInitialization[0][1]
    = InitializeRightFghHandler::getHandler(
    &simulation->_grid,
    &simulation->_parallelTopology,
    new RightMovingWallFghAction(parameters));

  simulation->_ghostCellsHandler._fghInitialization[1][0]
    = InitializeBottomFghHandler::getHandler(
    &simulation->_grid,
    &simulation->_parallelTopology,
    new BottomMovingWallFghAction(parameters));

  simulation->_ghostCellsHandler._fghInitialization[1][1]
    = InitializeTopFghHandler::getHandler(
    &simulation->_grid,
    &simulation->_parallelTopology,
    new TopMovingWallFghAction(parameters));

  if (D == 3) {
    simulation->_ghostCellsHandler._fghInitialization[2][0]
      = InitializeBackFghHandler::getHandler(
      &simulation->_grid,
      &simulation->_parallelTopology,
      new BackMovingWallFghAction(parameters));

    simulation->_ghostCellsHandler._fghInitialization[2][1]
      = InitializeFrontFghHandler::getHandler(
      &simulation->_grid,
      &simulation->_parallelTopology,
      new FrontMovingWallFghAction(parameters));
  }

  simulation->_ghostCellsHandler._velocityInitialization[0][0]
    = InitializeLeftVelocityHandler::getHandler(
    &simulation->_grid,
    &simulation->_parallelTopology,
    new LeftMovingWallVelocityAction(parameters,
                                     simulation->_maxVelocity));

  simulation->_ghostCellsHandler._velocityInitialization[0][1]
    = InitializeRightVelocityHandler::getHandler(
    &simulation->_grid,
    &simulation->_parallelTopology,
    new RightMovingWallVelocityAction(parameters,
                                      simulation->_maxVelocity));

  simulation->_ghostCellsHandler._velocityInitialization[1][0]
    = InitializeBottomVelocityHandler::getHandler(
    &simulation->_grid,
    &simulation->_parallelTopology,
    new BottomMovingWallVelocityAction(parameters,
                                       simulation->_maxVelocity));

  simulation->_ghostCellsHandler._velocityInitialization[1][1]
    = InitializeTopVelocityHandler::getHandler(
    &simulation->_grid,
    &simulation->_parallelTopology,
    new TopMovingWallVelocityAction(parameters,
                                    simulation->_maxVelocity));

  if (D == 3) {
    simulation->_ghostCellsHandler._velocityInitialization[2][0]
      = InitializeBackVelocityHandler::getHandler(
      &simulation->_grid,
      &simulation->_parallelTopology,
      new BackMovingWallVelocityAction(parameters,
                                       simulation->_maxVelocity));

    simulation->_ghostCellsHandler._velocityInitialization[2][1]
      = InitializeFrontVelocityHandler::getHandler(
      &simulation->_grid,
      &simulation->_parallelTopology,
      new FrontMovingWallVelocityAction(parameters,
                                        simulation->_maxVelocity));
  }

  return simulation;
}
}
}
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridFloat2D(Parameters& parameters) {
  return Private::createUniformGridFromTemplate<float, 2>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridDouble2D(Parameters& parameters) {
  return Private::createUniformGridFromTemplate<double, 2>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridFloat3D(Parameters& parameters) {
  return Private::createUniformGridFromTemplate<float, 3>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridDouble3D(Parameters& parameters) {
  return Private::createUniformGridFromTemplate<double, 3>(parameters);
}
