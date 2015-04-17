#pragma once

#include <Uni/ExecutionControl/exception>
#include <Uni/Logging/macros>
#include <Uni/math>

#include <precice/Constants.hpp>
#include <precice/SolverInterface.hpp>

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
template <unsigned TDimensions>
inline unsigned
get_distance(int const& position) {
  return (Uni::sign(position) * position) >> 2 * TDimensions;
}

template <unsigned TDimensions>
inline void
set_distance(int& position, unsigned const& distance) {
  position = Uni::sign(position) * (distance << 2 * TDimensions);
}

template <unsigned TDimensions>
inline unsigned
get_neighbor_bits(int const& position) {
  return (~(0 << (2 * TDimensions - 1))) & position;
}

template <unsigned TDimensions>
inline void
set_neighbor_bits(int& position, unsigned const& bits) {
  position = position | bits;
}

inline int
convert_precice_position(int const& position) {
  namespace pc = precice::constants;

  if (position == pc::positionOutsideOfGeometry()) {
    return 1;
  } else if (position == pc::positionOnGeometry()) {
    return -1;
  } else if (position == pc::positionInsideOfGeometry()) {
    return -1;
  }

  throwException("Failed to locate the position '{1}'",
                 position);

  return 0;
}

template <typename TCellAccessor>
inline void
set_cell_neighbors_along_geometry_interface(
  TCellAccessor const&      accessor,
  precice::SolverInterface* precice_interface,
  std::set<int> const&      body_meshes,
  unsigned const&           max_distance,
  unsigned const&           dimension) {
  auto const position = accessor.positionInRespectToGeometry()(dimension);

  auto const velocity_position = accessor.velocityPosition(dimension);

  unsigned distance = max_distance + 1;
  unsigned bits     = 0;

  for (unsigned d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (unsigned d2 = 0; d2 < 2; ++d2) {
      auto current_velocity_position = velocity_position;

      for (unsigned currentDistance = 1;
           currentDistance <= max_distance;
           ++currentDistance) {
        if (currentDistance > distance) {
          break;
        }

        int direction = -1;

        if (d2 == 1) {
          direction = +1;
        }

        // auto index  = accessor.index();
        // index(d) += direction * currentDistance;

        // if ((index(d)
        // < accessor.memory()->grid()->innerGrid.leftIndent() (d))
        // || (index(d)
        // >= accessor.memory()->grid()->innerGrid.innerLimit() (d))) {
        // continue;
        // }

        if (d == dimension) {
          current_velocity_position(dimension)
            += accessor.width(dimension, direction * currentDistance,
                              dimension);
        } else {
          current_velocity_position(d)
            += (accessor.width(d, direction * currentDistance, d)
                + accessor.width(d)) / 2.0;
        }

        auto const position1 = convert_precice_position(
          precice_interface->inquirePosition(
            current_velocity_position.template cast<double>().data(),
            body_meshes));

        if (Uni::sign(position) != Uni::sign(position1)) {
          if (currentDistance < distance) {
            bits     = 0;
            distance = currentDistance;
          }
          bits |= ((1 << d) << d2);
          break;
        }
      }
    }
  }

  set_distance<TCellAccessor::Dimensions>(
    accessor.positionInRespectToGeometry()(dimension),
    distance);
  set_neighbor_bits<TCellAccessor::Dimensions>(
    accessor.positionInRespectToGeometry()(dimension),
    bits);
}

template <typename TCellAccessor>
inline bool
validate_layer_number(TCellAccessor const& accessor,
                      unsigned const&      outer_layer_size,
                      unsigned const&      inner_layer_size,
                      unsigned const&      dimension) {
  auto const position = accessor.positionInRespectToGeometry()(dimension);

  auto const distance= get_distance<TCellAccessor::Dimensions>(position);

  if (position > 0) {
    if (outer_layer_size >= distance) {
      return true;
    }
  } else {
    if (inner_layer_size >= distance) {
      return true;
    }
  }
  return false;
}
}
}
}
