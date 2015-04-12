#include "SfsfdMemory.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template class SfsfdMemory
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 0, double, 2 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 1, double, 2 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 0, double, 2 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 1, double, 2 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 0, double, 3 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 1, double, 3 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 0, double, 3 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 1, double, 3 >>;
}
}
