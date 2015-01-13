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
template <typename TScalar, int TD>
SimulationFactory::Simulation*
createUniformGridFromTemplate(FluidSimulation::Configuration* parameters) {
  SimulationBuilder<TScalar, TD> builder(parameters);

  builder.setLeftAsParabolicInput();
  builder.setRightAsOutput();
  builder.setBottomAsMoving();
  builder.setTopAsMoving();

  if (TD == 3) {
    builder.setBackAsMoving();
    builder.setFrontAsMoving();
  }

  auto simulation = builder.simulation();

  return simulation;
}
}
}
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridFloat2D(FluidSimulation::Configuration* parameters) {
  return Private::createUniformGridFromTemplate<double, 2>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridDouble2D(FluidSimulation::Configuration* parameters) {
  return Private::createUniformGridFromTemplate<double, 2>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridFloat3D(FluidSimulation::Configuration* parameters) {
  return Private::createUniformGridFromTemplate<double, 2>(parameters);
}

SimulationFactory::Simulation*
SimulationFactory::
createUniformGridDouble3D(FluidSimulation::Configuration* parameters) {
  return Private::createUniformGridFromTemplate<double, 2>(parameters);
}
