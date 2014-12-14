#ifndef FsiSimulation_FluidSimulation_GhostLayer_Private_utilities_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_Private_utilities_hpp

#include <array>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <typename TFunctor, int TD>
using FunctorStack = std::array<std::array<TFunctor, 2>, TD>;
}
}
}

#endif
