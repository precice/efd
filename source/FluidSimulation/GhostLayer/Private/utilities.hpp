#ifndef FsiSimulation_Solvers_Ghost_handlerutilities_hpp
#define FsiSimulation_Solvers_Ghost_handlerutilities_hpp

#include <array>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <typename TFunctor, int D>
using FunctorStack = std::array<std::array<TFunctor, 2>, D>;
}
}
}

#endif
