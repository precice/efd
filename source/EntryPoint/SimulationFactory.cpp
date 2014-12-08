#include "SimulationFactory.hpp"

#include "Cell.hpp"
#include "GridGeometry.hpp"
#include "MyTemplateSimulation.hpp"
#include "Solvers/GhostInitializationActions.hpp"
#include "StructuredMemory/Accessor.hpp"
#include "StructuredMemory/Memory.hpp"

#include "MySimulationBuilder.hpp"

#include <Uni/StructuredGrid/Basic/MultiIndex>
//

using FsiSimulation::EntryPoint::SimulationFactory;

namespace FsiSimulation {
namespace EntryPoint {
namespace Private {
template <typename Scalar, int D>
SimulationFactory::Simulation*
createUniformGridFromTemplate(Parameters& parameters) {
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

  MySimulationBuilder<Scalar, D> builder(parameters);

  builder.setLeftWallMoving();
  builder.setRightWallMoving();
  builder.setBottomWallMoving();
  builder.setTopWallMoving();
  builder.setBackWallMoving();
  builder.setFrontWallMoving();

  auto simulation = builder.simulation();

  return simulation;
}
}
}
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridFloat2D(Parameters& parameters) {
  return Private::createUniformGridFromTemplate<double, 3>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridDouble2D(Parameters& parameters) {
  return Private::createUniformGridFromTemplate<double, 3>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridFloat3D(Parameters& parameters) {
  return Private::createUniformGridFromTemplate<double, 3>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridDouble3D(Parameters& parameters) {
  return Private::createUniformGridFromTemplate<double, 3>(parameters);
}
