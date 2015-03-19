#ifndef FsiSimulation_FluidSimulation_GhostLayer_Private_utilities_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_Private_utilities_hpp

#include "FluidSimulation/Configuration.hpp"

#include <Uni/Logging/macros>

#include <array>
#include <cmath>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <typename TFunctor, int TD>
using FunctorStack = std::array<std::array<TFunctor, 2>, TD>;

template <typename TCellAccessor,
          typename TVelocity>
inline
typename TCellAccessor::ScalarType
computeParabolicInputVelocity(
  TCellAccessor const& accessor,
  TVelocity const&     initialVelocity,
  int const&           d) {
  typedef TCellAccessor                           CellAccessorType;
  typedef typename CellAccessorType::VectorDiType VectorDiType;
  typedef typename CellAccessorType::ScalarType   ScalarType;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

  ScalarType newVelocity;

  newVelocity = std::pow(4, (Dimensions - 1))
                * initialVelocity(d);

  for (int i = 0; i < Dimensions; ++i) {
    if (i != d) {
      newVelocity *= accessor.velocityPosition(d, i)
                     * (accessor.memory()->gridGeometry()->size() (i)
                        - accessor.velocityPosition(d, i))
                     / (accessor.memory()->gridGeometry()->size() (i)
                        * accessor.memory()->gridGeometry()->size() (i));
    }
  }

  return newVelocity;
}
}
}
}

#endif
