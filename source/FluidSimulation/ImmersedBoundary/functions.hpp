#pragma once

#include <Uni/Logging/macros>

#include <precice/Constants.hpp>
#include <precice/SolverInterface.hpp>

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
inline int
is_outside(int const& position) {
  namespace pc = precice::constants;

  if (position == pc::positionOutsideOfGeometry()) {
    return 1;
  } else if (position == pc::positionOnGeometry()) {
    return 0;
  } else if (position == pc::positionInsideOfGeometry()) {
    return 0;
  }

  return -1;
}

inline int
is_the_same_position(int const& position0,
                     int const& position1) {
  auto const& result0 = is_outside(position0);
  auto const& result1 = is_outside(position1);

  if ((result0 == -1)
      || (result1 == -1)) {
    return -1;
  }

  if (result0 == result1) {
    return 1;
  } else {
    return 0;
  }
}

template <typename TCellAccessor>
inline int
compute_cell_layer_along_geometry_interface(TCellAccessor const& accessor,
                                            int const&           distance) {
  int position = accessor.positionInRespectToGeometry();

  bool result = true;

  for (int currentDistance = 1;
       currentDistance <= distance;
       ++currentDistance) {
    for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        auto index = accessor.index();

        if (d2 == 0) {
          index(d) -= currentDistance;
        } else {
          index(d) += currentDistance;
        }

        if ((index(d)
             < accessor.memory()->grid()->innerGrid.leftIndent() (d))
            || (index(d)
                >= accessor.memory()->grid()->innerGrid.innerLimit() (d))) {
          continue;
        }

        auto const position1
          = accessor.absolutePositionInRespectToGeometry(index);

        if (is_the_same_position(position, position1) == 0) {
          result = false;

          return currentDistance;
        } else {
          result = true;
        }
      }
    }
  }

  return -1;
}

template <typename TCellAccessor>
inline bool
validate_layer_number(TCellAccessor const& accessor,
                      int const&           distance,
                      int const&           outside_width,
                      int const&           inside_width) {
  ((void)outside_width);
  ((void)inside_width);

  Eigen::Vector4i pattern({ 1, 1, 1, 1 });

  if (distance != -1) {
    for (int i = 0; i < 4; ++i) {
      auto const& position = pattern(i);

      if (!((i == distance)
            && (position == 1))) {
        continue;
      }

      if (i < 2) {
        if (is_outside(accessor.positionInRespectToGeometry()) == 1) {
          return true;
        }
      } else {
        if (is_outside(accessor.positionInRespectToGeometry()) == 0) {
          return true;
        }
      }

      return false;
    }
  }

  return false;
}
}
}
}
