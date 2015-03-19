#pragma once

#include "FluidSimulation/Configuration.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
class Solver;
}
namespace EntryPoint {
class SolverFactory {
public:
  typedef FluidSimulation::Solver Solver;

public:
  static Solver*
  createSimpleFdDouble2D(FluidSimulation::Configuration*);

  static Solver*
  createSimpleFdDouble3D(FluidSimulation::Configuration*);

  static Solver*
  createFractionalStepFdDouble2D(FluidSimulation::Configuration*);

  static Solver*
  createFractionalStepDouble3D(FluidSimulation::Configuration*);
};
}
}
