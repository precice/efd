#include "SimulationFactory.hpp"

#include "Cell.hpp"
#include "GridGeometry.hpp"
#include "MyTemplateSimulation.hpp"
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
SimulationFactory::Simulation*
createUniformGridFromTemplate(Parameters const& parameters) {
  typedef SpecializedSimulationT<Scalar, D>
    SpecializedSimulation;
  typedef typename SpecializedSimulation::SpecializedGrid
    Grid;
  typedef typename SpecializedSimulation::SpecializedParallelTopology::VectorDi
    VectorDi;
  typedef typename SpecializedSimulation::GridGeometry::VectorDs
    VectorDs;

  auto simulation = new SpecializedSimulation();

  VectorDi processorSize;
  VectorDi globalCellSize;
  VectorDs size;
  processorSize(0)  = parameters.parallel.numProcessors[0];
  globalCellSize(0) = parameters.geometry.sizeX;
  size(0)           = parameters.geometry.lengthX;
  processorSize(1)  = parameters.parallel.numProcessors[1];
  globalCellSize(1) = parameters.geometry.sizeY;
  size(1)           = parameters.geometry.lengthY;

  if (D == 3) {
    processorSize(2)  = parameters.parallel.numProcessors[2];
    globalCellSize(2) = parameters.geometry.sizeZ;
    size(2)           = parameters.geometry.lengthZ;
  }

  simulation->_parallelTopology.initialize(parameters.parallel.rank,
                                           processorSize,
                                           globalCellSize);

  simulation->_parameters.re()    = parameters.flow.Re;
  simulation->_parameters.gamma() = parameters.solver.gamma;
  simulation->_parameters.g(0)    = parameters.environment.gx;
  simulation->_parameters.g(1)    = parameters.environment.gy;

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

  auto mpiVelocityExchangeStack =
    MpiVelocityExchangeStack::create(
      &simulation->_grid,
      &simulation->_parallelTopology);

  for (int d = 0; d < D; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      if (simulation->_parallelTopology.neighbors[d][d2] >= 0) {
        simulation->_ghostCellsHandler._mpiVelocityExchangeStack[d][d2] =
          mpiVelocityExchangeStack[d][d2];
      }
    }
  }

  simulation->_gridGeometry.initialize(size,
                                       simulation->_parallelTopology.
                                       globalSize,
                                       simulation->
                                       _parallelTopology
                                       .corner);

  return simulation;
}
}
}
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridFloat2D(Parameters const& parameters) {
  return Private::createUniformGridFromTemplate<float, 2>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridDouble2D(Parameters const& parameters) {
  return Private::createUniformGridFromTemplate<double, 2>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridFloat3D(Parameters const& parameters) {
  return Private::createUniformGridFromTemplate<float, 3>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridDouble3D(Parameters const& parameters) {
  return Private::createUniformGridFromTemplate<double, 3>(parameters);
}
