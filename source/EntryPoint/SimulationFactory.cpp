#include "SimulationFactory.hpp"

#include "FluidSimulation/Cell.hpp"
#include "FluidSimulation/FdSimulation.hpp"
#include "FluidSimulation/GridGeometry.hpp"

#include "SimulationBuilder.hpp"
//

using FsiSimulation::EntryPoint::SimulationFactory;

namespace FsiSimulation {
namespace EntryPoint {
namespace Private {
template <typename TScalar, int TD, int TSolverType = 0>
SimulationFactory::Simulation*
createUniformGridFromTemplate(FluidSimulation::Configuration* configuration) {
  SimulationBuilder<TScalar, TD, TSolverType> builder(configuration);

  auto simulation = builder.simulation();

  return simulation;
}
}
}
}

SimulationFactory::Simulation*
SimulationFactory::
createSimpleFdDouble2D(FluidSimulation::Configuration* configuration) {
  return Private::createUniformGridFromTemplate<double, 2, 1>(configuration);
}

SimulationFactory::Simulation*
SimulationFactory::
createSimpleFdDouble3D(FluidSimulation::Configuration* configuration) {
  return Private::createUniformGridFromTemplate<double, 2, 1>(configuration);
}

SimulationFactory::Simulation*
SimulationFactory::
createFractionalStepFdDouble2D(FluidSimulation::Configuration* configuration) {
  return Private::createUniformGridFromTemplate<double, 2, 1>(configuration);
}

SimulationFactory::Simulation*
SimulationFactory::
createFractionalStepDouble3D(FluidSimulation::Configuration* configuration) {
  return Private::createUniformGridFromTemplate<double, 2, 1>(configuration);
}
