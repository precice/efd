#pragma once

#include <Uni/Logging/macros>

#include <array>
#include <cmath>

namespace Fluid {
namespace Simulation {
namespace GhostLayer {
template <typename TFunctor, int TD>
using BasicFunctorStack = std::array<std::array<TFunctor, 2>, TD>;

template <typename TScalar>
inline TScalar
compute_parabolic_input(TScalar const& position,
                        TScalar const& width) {
  using ScalarType = TScalar;

  ScalarType result = 4;

  result *= position * (width - position) / (width * width);

  return result;
}

template <typename TVectorDs>
inline void
compute_parabolic_input(
  TVectorDs const&               position,
  TVectorDs const&               width,
  int const&                     dimension,
  typename TVectorDs::Scalar& value) {
  using VectorDsType = TVectorDs;

  enum {
    Dimensions = VectorDsType::RowsAtCompileTime
  };

  for (int i = 0; i < Dimensions; ++i) {
    if (i != dimension) {
      value *= compute_parabolic_input(position(i), width(i));
    }
  }
}

template <typename TCellAccessor>
inline void
compute_parabolic_input_velocity(TCellAccessor const&                accessor,
                                 int const&                          dimension,
                                 typename TCellAccessor::ScalarType& velocity) {

  compute_parabolic_input(accessor.velocityPosition(dimension),
                          accessor.memory()->gridGeometry()->size(),
                          dimension,
                          velocity);
}
}
}
}
