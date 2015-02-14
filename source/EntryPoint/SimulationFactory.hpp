#ifndef FsiSimulation_EntryPoint_SimulationFactory
#define FsiSimulation_EntryPoint_SimulationFactory

#include "FluidSimulation/Configuration.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
class Simulation;
}
namespace EntryPoint {
class SimulationFactory {
public:
  typedef FluidSimulation::Simulation Simulation;

public:
  static Simulation*
  createSimpleFdDouble2D(FluidSimulation::Configuration*);

  static Simulation*
  createSimpleFdDouble3D(FluidSimulation::Configuration*);

  static Simulation*
  createFractionalStepFdDouble2D(FluidSimulation::Configuration*);

  static Simulation*
  createFractionalStepDouble3D(FluidSimulation::Configuration*);
};
}
}
#endif
