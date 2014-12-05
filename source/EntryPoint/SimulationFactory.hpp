#ifndef FsiSimulation_EntryPoint_SimulationFactory
#define FsiSimulation_EntryPoint_SimulationFactory

#include "Parameters.h"

namespace FsiSimulation {
class MySimulation;
namespace EntryPoint {
class SimulationFactory {
public:
  typedef FsiSimulation::MySimulation Simulation;

public:
  static Simulation*
  createUniformGridFloat2D(Parameters& parameters);

  static Simulation*
  createUniformGridDouble2D(Parameters& parameters);

  static Simulation*
  createUniformGridFloat3D(Parameters& parameters);

  static Simulation*
  createUniformGridDouble3D(Parameters& parameters);
};
}
}
#endif
