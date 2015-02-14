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
typename TCellAccessor::CellType::Velocity::Scalar
computeParabolicInputVelocity(
  TCellAccessor const& accessor,
  TVelocity const&     initialVelocity,
  int const&           d) {
  typedef TCellAccessor                           CellAccessorType;
  typedef typename CellAccessorType::CellType     CellType;
  typedef typename CellAccessorType::VectorDiType VectorDiType;
  typedef typename CellType::Scalar               Scalar;

  enum {
    Dimensions = CellType::Dimensions
  };

  Scalar newVelocity;

  newVelocity = std::pow(4, (Dimensions - 1))
                * initialVelocity(d);

  for (int i = 0; i < Dimensions; ++i) {
    if (i != d) {
      newVelocity *= accessor.currentVelocityPosition(d)(i)
                     * (accessor.gridGeometry()->size() (i)
                        - accessor.currentVelocityPosition(d)(i))
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
