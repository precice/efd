#pragma once

#include "BodyForce/functions.hpp"

#include "Simulation/Grid.hpp"
#include "Simulation/IfsfdCellAccessor.hpp"
#include "Simulation/IfsfdMemory.hpp"
#include "Simulation/Private/mpigenerics.hpp"
#include "Simulation/Reporter.hpp"
#include "Simulation/SfsfdCellAccessor.hpp"
#include "Simulation/SfsfdMemory.hpp"

#include "functions.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/ExecutionControl/exception>

#include <array>
#include <map>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
namespace Private {
template <typename TVector>
struct BodyInterfaceCell {
  using Vector = TVector;
  enum {
    Dimensions = Vector::RowsAtCompileTime
  };

  BodyInterfaceCell(int const& vertexId_)
    : vertexId(vertexId_) {}

  int    vertexId;
  Vector data;
};
}
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

  using BodyInterfaceCellType = Private::BodyInterfaceCell<VectorDsType>;

  using VertexIdMap
          = std::map < int, std::unique_ptr < BodyInterfaceCellType >>;

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

  typename VertexIdMap::iterator
  begin() {
    return _vertexIds[0].begin();
  }

  typename VertexIdMap::iterator
  end() {
    return _vertexIds[0].end();
  }

  void
  resetIntermediateData() {
    for (auto& cell : this->_vertexIds[0]) {
      cell.second->data = VectorDsType::Zero();
    }
  }

  void
  initialize(precice::SolverInterface* preciceInterface,
             MemoryType const*         memory) {
    ((void)preciceInterface);
    ((void)memory);
  }

  void
  computePositionInRespectToGeometry(CellAccessorType const& accessor) const {
    ((void)accessor);
  }

  void
  createFluidMeshVertex(CellAccessorType const& accessor) {
    ((void)accessor);
  }

  std::pair<typename VertexIdMap::iterator, bool>
  doesVertexExist(CellAccessorType const& accessor) {
    ((void)accessor);

    return std::make_pair(end(), false);
  }

  void
  setFluidVelocity(typename VertexIdMap::iterator const& it,
                   VectorDsType const&                   velocity) {
    ((void)it);
    ((void)velocity);
  }

  void
  writeFluidVelocities() {}

  void
  mapData(ScalarType const& time_step_size) {
    ((void)time_step_size);
  }

  void
  readFluidForces() {}

  VectorDsType
  getFluidForce(typename VertexIdMap::iterator const& it) {
    ((void)it);

    return VectorDsType::Zero();
  }

  VectorDsType
  computeBodyForceAt(unsigned, VectorDsType const&, MemoryType*) {
    return VectorDsType::Zero();
  }

  void
  computeBodyForce(MemoryType const*, Reporter*) {}

protected:
  std::array<VertexIdMap, 1> _vertexIds;
};

template <typename TSolverTraits>
class Controller : public BasicController<TSolverTraits> {
public:
  using SolverTraitsType = TSolverTraits;

  using BaseType = BasicController<SolverTraitsType>;

  using BodyInterfaceCellType = typename BaseType::BodyInterfaceCellType;

  using VertexIdMap = typename BaseType::VertexIdMap;

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

public:
  Controller() :
    BaseType(),
    _maxLayerSize(2),
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

  typename VertexIdMap::iterator
  begin() {
    return this->_vertexIds[0].begin();
  }

  typename VertexIdMap::iterator
  end() {
    return this->_vertexIds[0].end();
  }

  void
  initialize(precice::SolverInterface* preciceInterface,
             MemoryType const*         memory) {
    // logInfo("Controller is being initialized");

    _preciceInterface = preciceInterface;
    _memory           = memory;

    _preciceInterface->initialize();
    _preciceInterface->initializeData();

    if (!_preciceInterface->hasMesh("FluidMesh")) {
      throwException("Precice configuration does not have 'FluidMesh'");
    }

    if (!_preciceInterface->hasMesh("BodyMesh")) {
      throwException("Precice configuration does not have 'BodyMesh'");
    }

    _fluidMeshId[0]
      = _preciceInterface->getMeshID("FluidMesh");
    _bodyMeshId
      = _preciceInterface->getMeshID("BodyMesh");

    if (!_preciceInterface->hasData("Velocities", _fluidMeshId[0])) {
      throwException("Precice configuration does not have 'Velocities' data"
                     " related to 'FluidMesh'");
    }

    if (!_preciceInterface->hasData("Forces", _fluidMeshId[0])) {
      throwException("Precice configuration does not have 'Forces' data"
                     " related to 'FluidMesh'");
    }

    _fluidMeshVelocitiesId[0]
      = _preciceInterface->getDataID("Velocities", _fluidMeshId[0]);
    _fluidMeshForcesId[0]
      = _preciceInterface->getDataID("Forces", _fluidMeshId[0]);
  }

  void
  computePositionInRespectToGeometry(CellAccessorType const& accessor) const {
    for (unsigned d = 0; d < Dimensions; ++d) {
      accessor.positionInRespectToGeometry(d)
        = convert_precice_position(
        _preciceInterface->inquirePosition(
          accessor.velocityPosition(d).template cast<double>().data(),
          std::set<int>({ _bodyMeshId })));
    }
    accessor.positionInRespectToGeometry(Dimensions)
      = convert_precice_position(
      _preciceInterface->inquirePosition(
        accessor.pressurePosition().template cast<double>().data(),
        std::set<int>({ _bodyMeshId })));
  }

  void
  createFluidMeshVertex(CellAccessorType const& accessor) {
    set_cell_neighbors_along_geometry_interface(
      accessor,
      _preciceInterface,
      std::set<int>({ _bodyMeshId }),
      _maxLayerSize);

    bool doAdd = validate_layer_number(accessor,
                                       outerLayerSize(),
                                       innerLayerSize());

    if (doAdd) {
      // logInfo("{1} | d = {2}", accessor.index().transpose(), distance);
      VectorDsType position = accessor.pressurePosition();

      auto vertexId
        = _preciceInterface->setMeshVertex(
        _fluidMeshId[0], position.data());

      this->_vertexIds[0].insert(
        std::make_pair(accessor.globalIndex(),
                       std::unique_ptr<BodyInterfaceCellType>(
                         new BodyInterfaceCellType(vertexId))));
    }
  }

  std::pair<typename VertexIdMap::iterator, bool>
  doesVertexExist(CellAccessorType const& accessor) {
    auto find_it = this->_vertexIds[0].find(accessor.globalIndex());

    if (find_it == this->_vertexIds[0].end()) {
      return std::make_pair(find_it, false);
    } else {
      return std::make_pair(find_it, true);
    }
  }

  void
  setFluidVelocity(typename VertexIdMap::iterator const& it,
                   VectorDsType const&                   velocity) {
    it->second->data = velocity;
  }

  void
  writeFluidVelocities() {
    for (auto const& cell : this->_vertexIds[0]) {
      VectorDsType resulting_velocity = cell.second->data;

      for (unsigned d = 0; d < Dimensions; ++d) {
        auto accessor = *_memory->grid()->innerGrid.begin();
        accessor.initialize(cell.first);
        auto neighbor_global_index = accessor.relativeGlobalIndex(d, -1);

        auto find_it = this->_vertexIds[0].find(neighbor_global_index);

        if (find_it == this->_vertexIds[0].end()) {
          // resulting_velocity(d) += 0.0;
        } else {
          resulting_velocity(d) += find_it->second->data(d);
        }
      }
      resulting_velocity /= 2.0;
      _preciceInterface->writeVectorData(
        _fluidMeshVelocitiesId[0],
        cell.second->vertexId,
        resulting_velocity.template cast<double>().data());
    }
  }

  void
  mapData(ScalarType const& time_step_size) {
    _preciceInterface->advance(time_step_size);
  }

  void
  readFluidForces() {
    this->resetIntermediateData();

    for (auto& cell : this->_vertexIds[0]) {
      Eigen::Matrix<double, Dimensions, 1> temp;
      _preciceInterface->readVectorData(
        _fluidMeshForcesId[0],
        cell.second->vertexId,
        temp.data());
      cell.second->data = temp.template cast<ScalarType>();
    }
  }

  VectorDsType
  getFluidForce(typename VertexIdMap::iterator const& it) {
    VectorDsType force = it->second->data;

    for (unsigned d = 0; d < Dimensions; ++d) {
      auto accessor = *_memory->grid()->innerGrid.begin();
      accessor.initialize(it->first);

      auto neighbor_global_index = accessor.relativeGlobalIndex(d, +1);

      auto find_it = this->_vertexIds[0].find(neighbor_global_index);

      if (find_it == this->_vertexIds[0].end()) {
        // force(d) += 0.0;
      } else {
        force(d) += find_it->second->data(d);
      }
    }
    force /= 2.0;

    return force / _memory->timeStepSize();
  }

  VectorDsType
  computeBodyForceAt(unsigned            global_index,
                     VectorDsType const& force,
                     MemoryType*         memory) {
    VectorDsType body_force = VectorDsType::Zero();

    auto accessor = *memory->grid()->begin();
    accessor.initialize(global_index);

    body_force += accessor.width().prod() * (-force);

    memory->addForceAt(global_index, body_force);

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

    FluidSimulation::Private::mpiAllReduce<ScalarType>(MPI_IN_PLACE,
                                                       total_force.data(),
                                                       Dimensions,
                                                       MPI_SUM,
                                                       PETSC_COMM_WORLD);

    FluidSimulation::Private::mpiAllReduce<ScalarType>(MPI_IN_PLACE,
                                                       total_force_turek.data(),
                                                       Dimensions,
                                                       MPI_SUM,
                                                       PETSC_COMM_WORLD);

    reporter->addAt(4, total_force);
    reporter->addAt(5, total_force_turek);
  }

protected:
  precice::SolverInterface* _preciceInterface;
  MemoryType const*         _memory;
  unsigned                  _maxLayerSize;
  unsigned                  _outerLayerSize;
  unsigned                  _innerLayerSize;

  std::array<int, Dimensions> _fluidMeshId;
  int                         _bodyMeshId;
  std::array<int, Dimensions> _fluidMeshVelocitiesId;
  std::array<int, Dimensions> _fluidMeshForcesId;
};
}
}
}
