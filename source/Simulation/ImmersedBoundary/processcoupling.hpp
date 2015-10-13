#pragma once

#include "Controller.hpp"
#include "functions.hpp"

#include "Simulation/Configuration.hpp"
#include "Simulation/functions.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <mpi.h>

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

template <typename TCellAccessor>
typename TCellAccessor::VectorDsType
compute_coupling_stress(
  TCellAccessor const&                      accessor,
  typename TCellAccessor::ScalarType const& diffusion_multiplier,
  typename TCellAccessor::ScalarType const& grad_pressure_multiplier) {
  using CellAccessor =  TCellAccessor;

  using Scalar = typename CellAccessor::ScalarType;

  using Vector = typename CellAccessor::VectorDsType;

  using Matrix = Eigen::Matrix<Scalar,
                               CellAccessor::Dimensions,
                               CellAccessor::Dimensions>;
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
        matrix(d, d2) = 0.0;
        logWarning("One of the neighboring cells in the dimension {1} "
                   "must be ouside the body; position {2}",
                   d2,
                   accessor.pressurePosition().transpose());
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

      Scalar width = 1.0;

      for (int d4 = 0; d4 < TCellAccessor::Dimensions; ++d4) {
        if (d4 != d) {
          width *= accessor.width(d4);
        }
      }

      matrix  = diffusion_multiplier * (matrix + matrix.transpose());
      matrix -= grad_pressure_multiplier * accessor.pressure() * Matrix::Identity();

      force += (matrix * normal) * width;
    }
  }

  return force;
}

template <typename TSolverTraits>
class CouplingController {
public:
  using SolverTraitsType = TSolverTraits;

  enum {Dimensions = SolverTraitsType::Dimensions};

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

public:
  CouplingController(
    typename FluidSimulation::Configuration const* configuration) {
    _couplingForcesName
      = configuration->get<std::string>("/Ib/Options/CouplingForcesName");
  }

  CouplingController(CouplingController const&) = delete;

  ~CouplingController() {}

  CouplingController&
  operator=(CouplingController const&) = delete;

  void
  initialize(MemoryType const*         memory,
             precice::SolverInterface* precice_interface) {
    _memory           = memory;
    _preciceInterface = precice_interface;

    if (!_preciceInterface->hasMesh("CouplingFluidMesh")) {
      throwException("Precice configuration does not have 'CouplingFluidMesh'");
    }
    _fluidMeshId = _preciceInterface->getMeshID("CouplingFluidMesh");

    if (!_preciceInterface->hasData(_couplingForcesName, _fluidMeshId)) {
      throwException("Precice configuration does not have '{1}' data"
                     " related to 'CouplingFluidMesh'",
                     _couplingForcesName);
    }

    _fluidMeshForcesId = _preciceInterface->getDataID(_couplingForcesName, _fluidMeshId);
  }

  void
  createMesh() {
    _preciceInterface->resetMesh(_fluidMeshId);
    _vertexIds.clear();

    std::vector<double> vertex_coords;

    unsigned index = 0;

    for (auto const& accessor : _memory->grid()->innerGrid) {
      if (!do_cell_force_computation(accessor)) {
        continue;
      }

      auto position = accessor.pressurePosition().template cast<double>();

      for (unsigned d = 0; d < Dimensions; ++d) {
        vertex_coords.push_back(position(d));
      }
      // _vertexIds.push_back(id++);
      ++index;
    }

    for (unsigned d = 0; d < Dimensions; ++d) {
      for (unsigned d2 = 0; d2 < 2; ++d2) {
        for (auto const& accessor : _memory->grid()->indentedBoundaries[d][d2]) {
          auto const positions = accessor.positionInRespectToGeometry();

          if (positions(Dimensions) != +1) {
            continue;
          }

          auto position = accessor.pressurePosition().template cast<double>();

          for (unsigned d3 = 0; d3 < Dimensions; ++d3) {
            vertex_coords.push_back(position(d3));
          }
          // _vertexIds.push_back(id++);
          ++index;
        }
      }
    }

    // unsigned id_offset = 0;

    // if ((_memory->parallelDistribution()->rank() - 1) >= 0) {
    // MPI_Status status;
    // MPI_Recv(&id_offset,
    // 1,
    // MPI_UNSIGNED,
    // (_memory->parallelDistribution()->rank() - 1),
    // 50,
    // _memory->parallelDistribution()->mpiCommunicator,
    // &status);
    // }

    // if ((_memory->parallelDistribution()->rank() + 1)
    // < _memory->parallelDistribution()->rankSize()) {
    // unsigned next_id_offset = id_offset + _vertexIds.size();

    // MPI_Send(&next_id_offset,
    // 1,
    // MPI_UNSIGNED,
    // (_memory->parallelDistribution()->rank() + 1),
    // 50,
    // _memory->parallelDistribution()->mpiCommunicator);
    // }

    // for (unsigned i = 0; i < _vertexIds.size(); ++i) {
    // _vertexIds[i] += id_offset;
    // }

    _vertexIds.resize(index);

    if (index > 0) {
      _preciceInterface->setMeshVertices(
        _fluidMeshId,
        index,
        vertex_coords.data(),
        _vertexIds.data());
    }

    // logInfo("Vertices {1} {2}", id, vertex_ids.size());
  }

  void
  sendStresses() {
    std::vector<double> stresses;
    stresses.resize(Dimensions * _vertexIds.size());

    unsigned index = 0;

    for (auto const& accessor : _memory->grid()->innerGrid) {
      if (!do_cell_force_computation(accessor)) {
        continue;
      }

      auto const force = compute_coupling_stress(
        accessor,
        _memory->parameters()->diffusionMultiplier(),
        _memory->parameters()->gradPressureMultiplier());

      for (unsigned d = 0; d < Dimensions; ++d) {
        stresses[index++] = static_cast<double>(force(d));
      }
      // logInfo("{1} {2}", id, temp_force.transpose());
    }

    for (unsigned d = 0; d < Dimensions; ++d) {
      for (unsigned d2 = 0; d2 < 2; ++d2) {
        for (auto const& accessor : _memory->grid()->indentedBoundaries[d][d2]) {
          auto const positions = accessor.positionInRespectToGeometry();

          if (positions(Dimensions) != +1) {
            continue;
          }

          // logInfo("{1}", accessor.pressurePosition());

          for (unsigned d3 = 0; d3 < Dimensions; ++d3) {
            stresses[index++] = 0.0;
          }
        }
      }
    }

    if (_vertexIds.size() > 0) {
      _preciceInterface->writeBlockVectorData(_fluidMeshForcesId,
                                              _vertexIds.size(),
                                              _vertexIds.data(),
                                              stresses.data());
    }
  }

private:
  MemoryType const*         _memory;
  precice::SolverInterface* _preciceInterface;
  std::string               _couplingForcesName;
  int                       _fluidMeshId;
  int                       _fluidMeshForcesId;

  std::vector<int> _vertexIds;
};
}
}
}
