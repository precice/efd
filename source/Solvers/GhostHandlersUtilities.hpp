#ifndef FsiSimulation_Solvers_Ghost_handlerutilities_hpp
#define FsiSimulation_Solvers_Ghost_handlerutilities_hpp

#include <array>

namespace FsiSimulation {
namespace Solvers {
namespace Ghost {
template <typename TFunctor, int D>
using FunctorStack = std::array<std::array<TFunctor, 2>, D>;
}
}
}

#endif
