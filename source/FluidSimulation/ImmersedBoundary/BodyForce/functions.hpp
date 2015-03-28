#pragma once

#include "FluidSimulation/ImmersedBoundary/functions.hpp"

#include "FluidSimulation/functions.hpp"

#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <precice/SolverInterface.hpp>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
namespace BodyForce {
template <typename TCellAccessor>
bool
computeCellForce(TCellAccessor const&                      accessor,
                 typename TCellAccessor::ScalarType const& re,
                 typename TCellAccessor::VectorDsType&     force) {
  typedef Eigen::Matrix<typename TCellAccessor::ScalarType,
                        TCellAccessor::Dimensions,
                        TCellAccessor::Dimensions> Matrix;

  typedef typename TCellAccessor::VectorDsType Vector;

  bool isCurrentOutside  = true;
  bool isNeighborsInside = false;

  if (is_outside(accessor.positionInRespectToGeometry()) == 1) {
    return false;
  }

  auto const layer_size
    = compute_cell_layer_along_geometry_interface(accessor, 2);

  if (layer_size != 2) {
    return false;
  }

  Matrix matrix;

  for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (int d2 = 0; d2 < TCellAccessor::Dimensions; ++d2) {
      matrix(d, d2) = dudx(
        accessor.velocity(d2, -1, d),
        accessor.velocity(d2, +1, d),
        accessor.width(d2, -1, d),
        accessor.width(d2, +1, d));
    }
  }

  matrix = (matrix + matrix.transpose()) / re ;

  matrix -= accessor.pressure() * Matrix::Identity();

  for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      int normal_direction = +1;
      int offset           = -2;

      if (d2 == 1) {
        normal_direction = -1;
        offset           = +2;
      }

      if (accessor.positionInRespectToGeometry(d, offset)
          != precice::constants::positionOutsideOfGeometry()) {
        Vector normal = Vector::Zero();

        normal(d) = normal_direction;

        Vector width = Vector::Zero();
        for (int d4 = 0; d4 < TCellAccessor::Dimensions; ++d4) {
          if (d4 != d) {
            width(d4) = normal(d) * accessor.width(d4);
          }
        }
        force = matrix * normal;
        // logInfo("Got it {1} {2} {3}", matrix, normal, force);
        break;
      }
    }
  }

  return true;
}
}
}
}
}
