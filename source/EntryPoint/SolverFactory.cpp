#include "SolverFactory.hpp"

#include "SolverBuilder.hpp"
//

using FsiSimulation::EntryPoint::SolverFactory;

namespace FsiSimulation {
namespace EntryPoint {
namespace Private {
template <typename TScalar, int TD, int TSolverType = 0>
SolverFactory::Solver*
createUniformGridFromTemplate(FluidSimulation::Configuration* configuration) {
  SolverBuilder<TScalar, TD, 0, 0> builder(configuration);

  auto simulation = builder.simulation();

  return simulation;
}
}
}
}

SolverFactory::Solver*
SolverFactory::
createSimpleFdDouble2D(FluidSimulation::Configuration* configuration) {
  return Private::createUniformGridFromTemplate<double, 2, 1>(configuration);
}

SolverFactory::Solver*
SolverFactory::
createSimpleFdDouble3D(FluidSimulation::Configuration* configuration) {
  return Private::createUniformGridFromTemplate<double, 2, 1>(configuration);
}

SolverFactory::Solver*
SolverFactory::
createFractionalStepFdDouble2D(FluidSimulation::Configuration* configuration) {
  return Private::createUniformGridFromTemplate<double, 2, 1>(configuration);
}

SolverFactory::Solver*
SolverFactory::
createFractionalStepDouble3D(FluidSimulation::Configuration* configuration) {
  return Private::createUniformGridFromTemplate<double, 2, 1>(configuration);
}
