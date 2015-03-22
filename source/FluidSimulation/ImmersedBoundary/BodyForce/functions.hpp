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
    = compute_cell_layer_along_geometry_interface(accessor, 1);

  if (layer_size == -1) {
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

  matrix = 1.0 / re * (matrix + matrix.transpose());

  matrix -= accessor.pressure() * Matrix::Identity();

  for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      int shift = -1;

      if (d2 == 1) {
        shift = 1;
      }

      if (accessor.positionInRespectToGeometry(d, shift)
          != precice::constants::positionOutsideOfGeometry()) {
        Vector normal = Vector::Zero();

        normal(d) = shift;

        for (int d4 = 0; d4 < TCellAccessor::Dimensions; ++d4) {
          if (d4 != d) {
            normal(d) *= accessor.width(d4);
          }
        }
        force += matrix * normal;
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
