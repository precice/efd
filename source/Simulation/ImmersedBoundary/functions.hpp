#pragma once

#include <Uni/ExecutionControl/exception>
#include <Uni/Logging/macros>
#include <Uni/math>

#include <precice/Constants.hpp>
#include <precice/SolverInterface.hpp>

#include <Eigen/Core>

namespace Fluid {
namespace Simulation {
namespace ImmersedBoundary {
template <unsigned TDimensions>
inline unsigned
get_distance(int const& location) {
  return Uni::sign(location) * location;
}

// template <unsigned TDimensions>
// inline unsigned
// get_distance(int const& location) {
// return (Uni::sign(location) * location) >> 2 * TDimensions;
// }

// template <unsigned TDimensions>
// inline unsigned
// get_neighbor_bits(int const& location) {
// return (~(0 << (2 * TDimensions - 1))) & location;
// }

// template <unsigned TDimensions>
// inline void
// set_position(int& location, unsigned const& distance, unsigned const& bits) {
// location = Uni::sign(location) * (distance << 2 * TDimensions);
// location = location | bits;
// }

inline int
convert_precice_position(int const& location) {
  namespace pc = precice::constants;

  if (location == pc::positionOutsideOfGeometry()) {
    return 1;
  } else if (location == pc::positionOnGeometry()) {
    return -1;
  } else if (location == pc::positionInsideOfGeometry()) {
    return -1;
  }

  throwException("Failed to find the location '{1}'", location);

  return 0;
}

template <typename TCellAccessor>
inline void
set_cell_neighbors_along_geometry_interface_in_dimension(
  TCellAccessor const&      accessor,
  precice::SolverInterface* precice_interface,
  std::set<int> const&      body_meshes,
  unsigned const&           max_distance,
  unsigned const&           dimension) {
  auto const location = accessor.positionInRespectToGeometry(dimension);

  typename TCellAccessor::VectorDsType position;

  if (dimension < TCellAccessor::Dimensions) {
    position = accessor.velocityPosition(dimension);
  } else {
    position = accessor.pressurePosition();
  }

  unsigned distance = max_distance + 1;
  // unsigned bits     = 0;

  for (unsigned d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (unsigned d2 = 0; d2 < 2; ++d2) {
      auto current_position = position;

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
          if (direction == -1) {
            current_position(d)
              -= accessor.width(d, -(currentDistance - 1), d);
          } else {
            current_position(d)
              += accessor.width(d, currentDistance, d);
          }
        } else {
          current_position(d)
            += direction
               * (accessor.width(d, direction * (currentDistance - 1), d)
                  + accessor.width(d, direction * currentDistance, d))
               / 2.0;
        }

        auto const position1 = convert_precice_position(
          precice_interface->inquirePosition(
            current_position.template cast<double>().data(),
            body_meshes));

        if (Uni::sign(location) != Uni::sign(position1)) {
          if (currentDistance < distance) {
            // bits     = 0;
            distance = currentDistance;
          }
          // bits |= ((1 << d) << d2);
          break;
        }
      }
    }
  }

  accessor.positionInRespectToGeometry(dimension)
    = Uni::sign(location) * distance;

  // set_position<TCellAccessor::Dimensions>(
  // accessor.positionInRespectToGeometry(dimension),
  // distance,
  // bits);
}

template <typename TCellAccessor>
inline void
set_cell_neighbors_along_geometry_interface(
  TCellAccessor const&      accessor,
  precice::SolverInterface* precice_interface,
  std::set<int> const&      body_meshes,
  unsigned const&           max_distance) {
  for (unsigned d = 0; d <= TCellAccessor::Dimensions; ++d) {
    set_cell_neighbors_along_geometry_interface_in_dimension(
      accessor,
      precice_interface,
      body_meshes,
      max_distance,
      d);
  }
}

template <typename TCellAccessor>
inline bool
validate_layer_number(TCellAccessor const& accessor,
                      unsigned const&      outer_layer_size,
                      unsigned const&      inner_layer_size) {
  auto const location
    = accessor.positionInRespectToGeometry(TCellAccessor::Dimensions);

  auto const distance = get_distance<TCellAccessor::Dimensions>(location);

  if (location > 0) {
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
