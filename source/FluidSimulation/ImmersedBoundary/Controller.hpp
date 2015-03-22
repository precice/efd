#pragma once

#include "functions.hpp"

#include <precice/SolverInterface.hpp>

#include <map>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
template <typename TSolverTraits>
class Controller {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridGeometryType = typename SolverTraitsType::GridGeometryType;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  using VertexIdMap = std::map<int, int>;

public:
  Controller() {}

  Controller(Controller const&) = delete;

  ~Controller() {}

  Controller&
  operator=(Controller const& other) = delete;

  typename VertexIdMap::const_iterator
  begin() const {
    return _vertexIds.begin();
  }

  typename VertexIdMap::const_iterator
  end() const {
    return _vertexIds.end();
  }

  void
  initialize(precice::SolverInterface* preciceInterface) {
    _preciceInterface = preciceInterface;

    _preciceInterface->initialize();
    _preciceInterface->initializeData();

    _fluidMeshId
      = _preciceInterface->getMeshID("FluidMesh");
    _fluidVertexSize
      = _preciceInterface->getMeshVertexSize(_fluidMeshId);
    _bodyMeshId
      = _preciceInterface->getMeshID("Body");
    _bodyVertexSize
      = _preciceInterface->getMeshVertexSize(_bodyMeshId);
    _fluidMeshVelocitiesId
      = _preciceInterface->getDataID("Velocities",
                                     _preciceInterface->getMeshID("FluidMesh"));
    _fluidMeshForcesId
      = _preciceInterface->getDataID("Forces",
                                     _preciceInterface->getMeshID("FluidMesh"));
    _bodyMeshVelocitiesId
      = _preciceInterface->getDataID("Velocities",
                                     _preciceInterface->getMeshID("Body"));
    _bodyMeshForcesId
      = _preciceInterface->getDataID("Forces",
                                     _preciceInterface->getMeshID("Body"));
  }

  void
  computePositionInRespectToGeometry(CellAccessorType const& accessor) const {
    accessor.positionInRespectToGeometry()
      = _preciceInterface->inquirePosition(
      accessor.pressurePosition().data(),
      std::set<int>({ _bodyMeshId }));
  }

  void
  createFluidMeshVertex(CellAccessorType const& accessor) {
    int distance
      = compute_cell_layer_along_geometry_interface(
      accessor, 2);

    bool doAdd = validate_layer_number(accessor, distance, 2, 2);

    if (doAdd) {
      VectorDsType position = accessor.pressurePosition();
      auto         vertexId
        = _preciceInterface->setMeshVertex(_fluidMeshId, position.data());

      _vertexIds.insert(std::make_pair(accessor.globalIndex(), vertexId));
    }
  }

  std::pair<typename VertexIdMap::const_iterator, bool>
  doesVertexExist(CellAccessorType const& accessor) const {
    auto find_it = _vertexIds.find(accessor.globalIndex());

    if (find_it == _vertexIds.end()) {
      return std::make_pair(find_it, false);
    } else {
      return std::make_pair(find_it, true);
    }
  }

  void
  writeFluidVelocity(typename VertexIdMap::const_iterator const& it,
                     VectorDsType const&                         velocity) {
    _preciceInterface->writeVectorData(_fluidMeshVelocitiesId,
                                       it->second,
                                       velocity.data());
  }

  void
  mapData(ScalarType const& time_step_size) {
    _preciceInterface->advance(time_step_size);
  }

  void
  readFluidForce(typename VertexIdMap::const_iterator const& it,
                 VectorDsType&                               force) {
    // double data[Dimensions];
    _preciceInterface->readVectorData(_fluidMeshForcesId, it->second,
                                      force.data());

    // force
    //   = Eigen::Matrix<double, Dimensions, 1>(data).template cast<ScalarType>();
  }

private:
  VertexIdMap _vertexIds;

  precice::SolverInterface* _preciceInterface;

  int _fluidMeshId;
  int _fluidVertexSize;
  int _bodyMeshId;
  int _bodyVertexSize;
  int _fluidMeshVelocitiesId;
  int _fluidMeshForcesId;
  int _bodyMeshVelocitiesId;
  int _bodyMeshForcesId;
};
}
}
}
