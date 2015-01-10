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

template <typename TCellAcessor,
          typename TVelocity>
inline
typename TCellAcessor::CellType::Velocity::Scalar
computeParabolicInputVelocity(
  TCellAcessor const& accessor,
  TVelocity const&    initialVelocity,
  int const&          d) {
  typedef typename TCellAcessor::CellType::Velocity::Scalar Scalar;

  Scalar newVelocity;

  newVelocity = std::pow(4, (accessor.gridGeometry()->size().size() - 1))
                * initialVelocity(d);

  for (int i = 0; i < accessor.gridGeometry()->size().size(); ++i) {
    if (i != d) {
      newVelocity *= accessor.currentVelocityPosition (d)(i)
                     * (accessor.gridGeometry()->size() (i)
                        - accessor.currentVelocityPosition (d)(i))
                     / (accessor.gridGeometry()->size() (i)
                        * accessor.gridGeometry()->size() (i));
    }
  }

  return newVelocity;
}
}
}
}

#endif
