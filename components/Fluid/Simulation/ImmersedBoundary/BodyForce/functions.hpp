#pragma once

#include "Simulation/ImmersedBoundary/functions.hpp"

#include "Simulation/functions.hpp"

#include <Uni/Logging/macros>

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
namespace BodyForce {
template <typename TCellAccessor>
inline Eigen::Matrix<int, 2* TCellAccessor::Dimensions, 1>
compute_neighbor_locations(TCellAccessor const& accessor) {
  using Vector = Eigen::Matrix<int, 2* TCellAccessor::Dimensions, 1>;
  Vector locations = Vector::Zero();

  for (unsigned d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      int direction_offset = -1;

      if (d2 == 1) {
        direction_offset = +1;
      }

      bool isOuside = true;

      for (unsigned d3 = 0; d3 <= TCellAccessor::Dimensions; ++d3) {
        if (accessor.positionInRespectToGeometry(d, direction_offset, d3) < 0) {
          isOuside = false;
          break;
        }
      }

      if (isOuside) {
        locations(2 * d + d2) = 1;
      }
    }
  }

  return locations;
}

template <typename TCellAccessor>
bool
do_cell_force_computation(TCellAccessor const& accessor) {
  auto const positions = accessor.positionInRespectToGeometry();

  bool doAccept = false;

  for (unsigned d = 0; d <= TCellAccessor::Dimensions; ++d) {
    if (positions(d) < 0) {
      return false;
    }

    if (get_distance<TCellAccessor::Dimensions>(positions(d)) == 1) {
      doAccept = true;
    }
  }

  return doAccept;
}

template <typename TCellAccessor>
bool
compute_cell_force(TCellAccessor const&                      accessor,
                   typename TCellAccessor::ScalarType const& re,
                   typename TCellAccessor::VectorDsType&     force) {
  using Scalar = typename TCellAccessor::ScalarType;

  using Vector = typename TCellAccessor::VectorDsType;

  using Matrix = Eigen::Matrix<Scalar,
                               TCellAccessor::Dimensions,
                               TCellAccessor::Dimensions>;

  if (!do_cell_force_computation(accessor)) {
    return false;
  }

  auto locations = compute_neighbor_locations(accessor);

  Matrix matrix;

  for (unsigned d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (unsigned d2 = 0; d2 < TCellAccessor::Dimensions; ++d2) {
      if (locations(2 * d2 + 0) == 1
          && locations(2 * d2 + 1) == 1) {
        matrix(d, d2)
          = (accessor.velocity(d2, +1, d) - accessor.velocity(d2, -1, d))
            / (accessor.width(d2)
               + (accessor.width(d2, +1, d2) + accessor.width(d2, -1, d2))
               / 2.0);
      } else if (locations(2 * d2 + 0) == 1) {
        matrix(d, d2)
          = (accessor.velocity(d) - accessor.velocity(d2, -1, d))
            / accessor.width(d2);
      } else if (locations(2 * d2 + 1) == 1) {
        matrix(d, d2)
          = (accessor.velocity(d2, +1, d) - accessor.velocity(d))
            / accessor.width(d2, +1, d2);
      } else {
        throwException("One of the neighboring cells in the dimension {1} "
                       "must be ouside the body", d2);
      }
    }
  }

  force = Vector::Zero();

  for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      if (locations(2 * d + d2) == 1) {
        continue;
      }

      int normal_direction = +1.0;

      if (d2 == 1) {
        normal_direction = -1.0;
      }

      Vector normal = Vector::Zero();

      normal(d) = normal_direction;

      Scalar width = 1.0;

      for (int d4 = 0; d4 < TCellAccessor::Dimensions; ++d4) {
        if (d4 != d) {
          width *= accessor.width(d4);
        }
      }

      matrix = (matrix + matrix.transpose()) / re;

      matrix -= accessor.pressure() * Matrix::Identity();

      normal(d) *= width;
      force     += matrix * normal;
    }
  }

  return true;
}

template <typename TCellAccessor>
bool
compute_cell_force_turek(TCellAccessor const&                      accessor,
                         typename TCellAccessor::ScalarType const& re,
                         typename TCellAccessor::VectorDsType&     force) {
  using Scalar = typename TCellAccessor::ScalarType;

  using Vector = typename TCellAccessor::VectorDsType;

  if (!do_cell_force_computation(accessor)) {
    return false;
  }

  auto locations = compute_neighbor_locations(accessor);

  force = Vector::Zero();

  for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      if (locations(2 * d + d2) == 1) {
        continue;
      }

      int normal_direction = +1.0;

      if (d2 == 1) {
        normal_direction = -1.0;
      }

      Vector normal = Vector::Zero();

      normal(d) = normal_direction;

      Scalar width = 1.0;

      for (int d3 = 0; d3 < TCellAccessor::Dimensions; ++d3) {
        if (d3 != d) {
          width *= accessor.width(d3);
        }
      }
      Vector tangent = normal;
      tangent(0) = normal(1);
      tangent(1) = normal(0);

      Vector grad_velocity;

      for (int d3 = 0; d3 < TCellAccessor::Dimensions; ++d3) {
        if (locations(2 * d3 + 0) == 1
            && locations(2 * d3 + 1) == 1) {
          grad_velocity(d3)
            = (accessor.velocity(d3, +1).dot(tangent)
               - accessor.velocity(d3, -1).dot(tangent))
              / (accessor.width(d3)
                 + (accessor.width(d3, +1, d3)
                    + accessor.width(d3, -1, d3)) / 2.0);
        } else if (locations(2 * d3 + 0) == 1) {
          grad_velocity(d3)
            = (accessor.velocity().dot(tangent)
               - accessor.velocity(d3, -1).dot(tangent))
              / accessor.width(d3);
        } else if (locations(2 * d3 + 1) == 1) {
          grad_velocity(d3)
            = (accessor.velocity(d3, +1).dot(tangent)
               - accessor.velocity().dot(tangent))
              / accessor.width(d3, +1, d3);
        } else {
          throwException("One of the neighboring cells in the dimension {1} "
                         "must be ouside the body", d3);
        }
      }

      Scalar dvdn = grad_velocity.dot(normal) / re;
      force(0)
        += (dvdn * normal(1) - accessor.pressure() * normal(0)) * width;
      force(1)
        += (dvdn * normal(0) - accessor.pressure() * normal(1)) * width;
    }
  }

  return true;
}
}
}
}
}
