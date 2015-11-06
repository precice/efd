#pragma once

#include <cmath>

namespace Fluid {
namespace Simulation {
template <typename TCellAccessor>
void
compute_max_velocity(TCellAccessor const&                  accessor,
                   typename TCellAccessor::VectorDsType& maxVelocity) {
  maxVelocity
    = accessor.velocity().cwiseAbs().cwiseQuotient(
    accessor.width()).cwiseMax(maxVelocity);
}

template <typename TVector>
void
compute_max_velocity(TVector const& velocity,
                     TVector const& width,
                     TVector&       maxVelocity) {
  maxVelocity
    = velocity.cwiseAbs().cwiseQuotient(width).cwiseMax(maxVelocity);
}

template <typename TVector, typename TScalar>
void
compute_max_velocity(TScalar const&  velocity,
                     TScalar const&  width,
                     unsigned const& dimension,
                     TVector&        maxVelocity) {
  maxVelocity(dimension)
    = std::max(std::abs(velocity) / width, maxVelocity(dimension));
}
}
}
