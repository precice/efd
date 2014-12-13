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
  createUniformGridFloat2D(FluidSimulation::Parameters& parameters);

  static Simulation*
  createUniformGridDouble2D(FluidSimulation::Parameters& parameters);

  static Simulation*
  createUniformGridFloat3D(FluidSimulation::Parameters& parameters);

  static Simulation*
  createUniformGridDouble3D(FluidSimulation::Parameters& parameters);
};
}
}
#endif
