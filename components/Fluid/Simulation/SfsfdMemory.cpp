#include "SfsfdMemory.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template class SfsfdMemory
  < SfsfdSolverTraits < 0, 0, double, 2 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < 0, 1, double, 2 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < 1, 0, double, 2 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < 1, 1, double, 2 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < 0, 0, double, 3 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < 0, 1, double, 3 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < 1, 0, double, 3 >>;
template class SfsfdMemory
  < SfsfdSolverTraits < 1, 1, double, 3 >>;
}
}
