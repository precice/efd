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
compute_cell_force(TCellAccessor const&                      accessor,
                   typename TCellAccessor::ScalarType const& re,
                   typename TCellAccessor::VectorDsType&     force) {
  typedef Eigen::Matrix<typename TCellAccessor::ScalarType,
                        TCellAccessor::Dimensions,
                        TCellAccessor::Dimensions> Matrix;

  using Scalar = typename TCellAccessor::ScalarType;
  using Vector = typename TCellAccessor::VectorDsType;

  bool isCurrentOutside  = true;
  bool isNeighborsInside = false;

  if (is_outside(accessor.positionInRespectToGeometry()) != 1) {
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
        accessor.width(d),
        accessor.width(d2, -1, d),
        accessor.width(d2, +1, d));
    }
  }

  matrix = (matrix + matrix.transpose()) / re;

  matrix -= accessor.pressure() * Matrix::Identity();

  force = Vector::Zero();

  for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      Scalar normal_direction = +1.0;
      int offset           = -2;

      if (d2 == 1) {
        normal_direction = -1.0;
        offset           = +2;
      }

      if (accessor.positionInRespectToGeometry(d, offset)
          != precice::constants::positionOutsideOfGeometry()) {
        Vector normal = Vector::Zero();

        normal(d) = normal_direction;

        Scalar width = 1.0;

        for (int d4 = 0; d4 < TCellAccessor::Dimensions; ++d4) {
          if (d4 != d) {
            width *= accessor.width(d4);
          }
        }
        normal(d) *= width;
        force     += matrix * normal;
      }
    }
  }

  return true;
}

template <typename TCellAccessor>
bool
compute_cell_force_turek(TCellAccessor const&                      accessor,
                         typename TCellAccessor::ScalarType const& re,
                         typename TCellAccessor::VectorDsType&     force) {
  typedef Eigen::Matrix<typename TCellAccessor::ScalarType,
                        TCellAccessor::Dimensions,
                        TCellAccessor::Dimensions> Matrix;

  using Scalar = typename TCellAccessor::ScalarType;
  using Vector = typename TCellAccessor::VectorDsType;

  bool isCurrentOutside  = true;
  bool isNeighborsInside = false;

  if (is_outside(accessor.positionInRespectToGeometry()) != 1) {
    return false;
  }

  auto const layer_size
    = compute_cell_layer_along_geometry_interface(accessor, 2);

  if (layer_size != 2) {
    return false;
  }

  force = Vector::Zero();

  for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      Scalar normal_direction = +1.0;
      int    offset           = -2;

      if (d2 == 1) {
        normal_direction = -1.0;
        offset           = +2;
      }

      if (accessor.positionInRespectToGeometry(d, offset)
          != precice::constants::positionOutsideOfGeometry()) {
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
          grad_velocity(d3) = dudx(
            accessor.velocity(d3, -1).dot(tangent),
            accessor.velocity(d3, +1).dot(tangent),
            accessor.width(d3),
            accessor.width(d3, -1, d3),
            accessor.width(d3, +1, d3));
        }

        Scalar dvdn = grad_velocity.dot(normal) / re;
        force(0)
          += (dvdn * normal(1) - accessor.pressure() * normal(0)) * width;
        force(1)
          += (dvdn * normal(0) - accessor.pressure() * normal(1)) * width;
      }
    }
  }

  return true;
}
}
}
}
}
