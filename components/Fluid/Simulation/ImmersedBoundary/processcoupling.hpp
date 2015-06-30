#pragma once

#include "functions.hpp"

#include "Simulation/functions.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/Logging/macros>

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
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

      bool is_outside = true;

      for (unsigned d3 = 0; d3 <= TCellAccessor::Dimensions; ++d3) {
        if (accessor.positionInRespectToGeometry(d, direction_offset, d3) < 0) {
          is_outside = false;
          break;
        }
      }

      if (is_outside) {
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

template <typename TSubgrid>
void
create_coupling_mesh(TSubgrid const&           subgrid,
                     precice::SolverInterface* preciceInterface) {
  using CellAccessor = typename TSubgrid::CellAccessor;

  using Scalar = typename CellAccessor::ScalarType;

  using Vector = typename CellAccessor::VectorDsType;

  if (!preciceInterface->hasMesh("FluidMesh")) {
    throwException("Precice configuration does not have 'FluidMesh'");
  }

  auto const fluidMeshId = preciceInterface->getMeshID("FluidMesh");

  preciceInterface->resetMesh(fluidMeshId);

  std::vector<int>    vertex_ids;
  std::vector<double> vertex_coords;

  unsigned id = 0;

  for (auto const& accessor : subgrid) {
    if (!do_cell_force_computation(accessor)) {
      continue;
    }

    auto position = accessor.pressurePosition().template cast<double>();

    for (unsigned d = 0; d < CellAccessor::Dimensions; ++d) {
      vertex_coords.push_back(position(d));
    }
    vertex_ids.push_back(id++);
  }
  preciceInterface->setMeshVertices(
    fluidMeshId,
    vertex_ids.size(),
    vertex_coords.data(),
    vertex_ids.data());
  // logInfo("Vertices {1} {2}", id, vertex_ids.size());
}

template <typename TSubgrid>
void
send_coupling_stresses(
  TSubgrid const&                                    subgrid,
  precice::SolverInterface*                          preciceInterface,
  typename TSubgrid::CellAccessor::ScalarType const& diffusion_multiplier,
  typename TSubgrid::CellAccessor::ScalarType const& grad_pressure_multiplier) {
  using CellAccessor = typename TSubgrid::CellAccessor;

  using Scalar = typename CellAccessor::ScalarType;

  using Vector = typename CellAccessor::VectorDsType;

  using Matrix = Eigen::Matrix<Scalar,
                               TSubgrid::Dimensions,
                               TSubgrid::Dimensions>;

  if (!preciceInterface->hasMesh("FluidMesh")) {
    throwException("Precice configuration does not have 'FluidMesh'");
  }
  auto const fluidMeshId = preciceInterface->getMeshID("FluidMesh");

  if (!preciceInterface->hasData("CouplingStresses", fluidMeshId)) {
    throwException("Precice configuration does not have 'CouplingStresses' data"
                   " related to 'FluidMesh'");
  }

  auto const fluidMeshStressesId = preciceInterface->getDataID("CouplingStresses",
                                                               fluidMeshId);

  unsigned id = 0;

  for (auto const& accessor : subgrid) {
    if (!do_cell_force_computation(accessor)) {
      continue;
    }

    auto locations = compute_neighbor_locations(accessor);

    Matrix matrix;

    for (unsigned d = 0; d < CellAccessor::Dimensions; ++d) {
      for (unsigned d2 = 0; d2 < CellAccessor::Dimensions; ++d2) {
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

    Vector force = Vector::Zero();

    for (int d = 0; d < CellAccessor::Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        if (locations(2 * d + d2) == 1) {
          continue;
        }

        int normal_direction = +1;

        if (d2 == 1) {
          normal_direction = -1;
        }

        Vector normal = Vector::Zero();

        normal(d) = normal_direction;

        matrix  = diffusion_multiplier * (matrix + matrix.transpose());
        matrix -= grad_pressure_multiplier * accessor.pressure() * Matrix::Identity();

        force += matrix * normal;
      }
    }
    auto const temp_force = force.template cast<double>().eval();
    preciceInterface->writeVectorData(fluidMeshStressesId,
                                      id++,
                                      temp_force.data());
    // logInfo("{1} {2}", id, temp_force.transpose());
  }
}
}
}
}
