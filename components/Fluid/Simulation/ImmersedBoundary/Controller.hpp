#pragma once

#include "BodyForce/functions.hpp"
#include "Simulation/Reporter.hpp"
#include "functions.hpp"

#include "Simulation/Private/mpigenerics.hpp"

#include <precice/SolverInterface.hpp>

#include <map>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
/**
 * TODO:
 *      1. Change vertexId from global index to unique id
 */
template <typename TSolverTraits>
class BasicController {
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
  BasicController() {}

  BasicController(BasicController const&) = delete;

  ~BasicController() {}

  BasicController&
  operator=(BasicController const& other) = delete;

  unsigned const&
  outerLayerSize() const {
    static unsigned value = 0;

    return value;
  }

  unsigned&
  outerLayerSize(unsigned const& i) {
    ((void)i);
    static unsigned value = 0;

    return value;
  }

  unsigned const&
  innerLayerSize() const {
    static unsigned value = 0;

    return value;
  }

  unsigned&
  innerLayerSize(unsigned const& i) {
    ((void)i);
    static unsigned value = 0;

    return value;
  }

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
    ((void)preciceInterface);
  }

  void
  computePositionInRespectToGeometry(CellAccessorType const& accessor) const {
    ((void)accessor);
  }

  void
  createFluidMeshVertex(CellAccessorType const& accessor) {
    ((void)accessor);
  }

  std::pair<typename VertexIdMap::const_iterator, bool>
  doesVertexExist(CellAccessorType const& accessor) const {
    ((void)accessor);
    return std::make_pair(end(), false);
  }

  void
  writeFluidVelocity(typename VertexIdMap::const_iterator const& it,
                     VectorDsType const&                         velocity) {
    ((void)it);
    ((void)velocity);
  }

  void
  mapData(ScalarType const& time_step_size) {
    ((void)time_step_size);
  }

  void
  readFluidForce(typename VertexIdMap::const_iterator const& it,
                 VectorDsType&                               force) {
    ((void)it);
    force = VectorDsType::Zero();
  }

  VectorDsType
  computeBodyForceAt(unsigned, VectorDsType const&, MemoryType*) {
    return VectorDsType::Zero();
  }

  void
  computeBodyForce(MemoryType const*, Reporter*) {}

protected:
  VertexIdMap _vertexIds;
};

template <typename TSolverTraits>
class Controller : public BasicController<TSolverTraits> {
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
  Controller() :
    _maxLayerSize(0),
    _outerLayerSize(0),
    _innerLayerSize(0) {}

  Controller(Controller const&) = delete;

  ~Controller() {}

  Controller&
  operator=(Controller const& other) = delete;

  unsigned const&
  outerLayerSize() const {
    return _outerLayerSize;
  }

  unsigned&
  outerLayerSize(unsigned const& i) {
    _maxLayerSize = std::max(_maxLayerSize, i);

    return _outerLayerSize = i;
  }

  unsigned const&
  innerLayerSize() const {
    return _innerLayerSize;
  }

  unsigned&
  innerLayerSize(unsigned const& i) {
    _maxLayerSize = std::max(_maxLayerSize, i);

    return _innerLayerSize = i;
  }

  typename VertexIdMap::const_iterator
  begin() const {
    return this->_vertexIds.begin();
  }

  typename VertexIdMap::const_iterator
  end() const {
    return this->_vertexIds.end();
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
      accessor,
      _maxLayerSize);

    bool doAdd = validate_layer_number(accessor,
                                       distance,
                                       outerLayerSize(),
                                       innerLayerSize());

    if (doAdd) {
      VectorDsType position = accessor.pressurePosition();
      auto         vertexId
        = _preciceInterface->setMeshVertex(_fluidMeshId, position.data());

      this->_vertexIds.insert(std::make_pair(accessor.globalIndex(), vertexId));
    }
  }

  std::pair<typename VertexIdMap::const_iterator, bool>
  doesVertexExist(CellAccessorType const& accessor) const {
    auto find_it = this->_vertexIds.find(accessor.globalIndex());

    if (find_it == this->_vertexIds.end()) {
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
    _preciceInterface->readVectorData(_fluidMeshForcesId,
                                      it->second,
                                      force.data());
  }

  VectorDsType
  computeBodyForceAt(unsigned            global_index,
                     VectorDsType const& force,
                     MemoryType*         memory) {
    VectorDsType body_force = VectorDsType::Zero();

    auto accessor = *memory->grid()->begin();
    accessor.initialize(global_index);

    body_force += accessor.width().prod() * (-force) / memory->timeStepSize();

    memory->setForceAt(global_index, body_force);

    return body_force;
  }

  void
  computeBodyForce(MemoryType const* memory,
                   Reporter*         reporter) {
    VectorDsType total_force       = VectorDsType::Zero();
    VectorDsType total_force_turek = VectorDsType::Zero();

    for (auto const& accessor : memory->grid()->innerGrid) {
      VectorDsType force       = VectorDsType::Zero();
      VectorDsType force_turek = VectorDsType::Zero();
      BodyForce::compute_cell_force(accessor,
                                    memory->parameters()->re(),
                                    force);
      BodyForce::compute_cell_force_turek(accessor,
                                          memory->parameters()->re(),
                                          force_turek);
      accessor.setBodyForce(force);
      total_force       += force;
      total_force_turek += force_turek;
    }

    Private::mpiAllReduce<ScalarType>(MPI_IN_PLACE,
                                      total_force.data(),
                                      Dimensions,
                                      MPI_SUM,
                                      PETSC_COMM_WORLD);

    Private::mpiAllReduce<ScalarType>(MPI_IN_PLACE,
                                      total_force_turek.data(),
                                      Dimensions,
                                      MPI_SUM,
                                      PETSC_COMM_WORLD);

    reporter->addAt(4, total_force);
    reporter->addAt(5, total_force_turek);
  }

private:
  precice::SolverInterface* _preciceInterface;
  unsigned                  _maxLayerSize;
  unsigned                  _outerLayerSize;
  unsigned                  _innerLayerSize;

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
