#ifndef FsiSimulation_EntryPoint_SimulationFactory
#define FsiSimulation_EntryPoint_SimulationFactory

#include "Cell.hpp"
#include "CellHeap.hpp"
#include "GridGeometry.hpp"
#include "MySimulation.hpp"

namespace FsiSimulation {
namespace EntryPoint {
class SimulationFactory {
public:
  typedef FsiSimulation::MySimulation Simulation;

public:
  static Simulation*
  createUniformGridFloat2D() {
    return createUniformGridFromTemplate<float, 2>();
  }

  static Simulation*
  createUniformGridDouble2D() {
    return createUniformGridFromTemplate<double, 2>();
  }

  static Simulation*
  createUniformGridFloat3D() {
    return createUniformGridFromTemplate<float, 3>();
  }

  static Simulation*
  createUniformGridDouble3D() {
    return createUniformGridFromTemplate<double, 3>();
  }

private:
  template <typename Scalar, int D>
  static
  Simulation*
  createUniformGridFromTemplate() {
    return new MyTemplateSimulation<UniformGridGeometry<Scalar, D>,
                                    CellHeap<Cell<Scalar, D>, D>,
                                    Scalar,
                                    D>();
  }
};
}
}
#endif
