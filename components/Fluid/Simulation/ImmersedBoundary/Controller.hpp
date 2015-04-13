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
  begin(unsigned const& dimension) const {
    return _vertexIds[dimension].begin();
  }

  typename VertexIdMap::const_iterator
  end(unsigned const& dimension) const {
    return _vertexIds[dimension].end();
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
  doesVertexExist(CellAccessorType const& accessor,
                  unsigned const&         dimension) const {
    ((void)accessor);

    return std::make_pair(end(dimension), false);
  }

  void
  writeFluidVelocity(typename VertexIdMap::const_iterator const& it,
                     ScalarType const&                           velocity,
                     unsigned const&                             dimension) {
    ((void)it);
    ((void)velocity);
    ((void)dimension);
  }

  void
  mapData(ScalarType const& time_step_size) {
    ((void)time_step_size);
  }

  void
  readFluidForce(typename VertexIdMap::const_iterator const& it,
                 ScalarType&                                 force,
                 unsigned const&                             dimension) {
    ((void)it);
    ((void)dimension);
    force = 0.0;
  }

  VectorDsType
  computeBodyForceAt(unsigned, ScalarType const&, MemoryType*,
                     unsigned const&) {
    return VectorDsType::Zero();
  }

  void
  computeBodyForce(MemoryType const*, Reporter*) {}

protected:
  std::array<VertexIdMap, Dimensions> _vertexIds;
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
  begin(unsigned const& dimension) const {
    return this->_vertexIds[dimension].begin();
  }

  typename VertexIdMap::const_iterator
  end(unsigned const& dimension) const {
    return this->_vertexIds[dimension].end();
  }

  void
  initialize(precice::SolverInterface* preciceInterface) {
    // logInfo("Controller is being initialized");

    _preciceInterface = preciceInterface;

    _preciceInterface->initialize();
    _preciceInterface->initializeData();

    if (!_preciceInterface->hasMesh("FluidMeshX")) {
      throwException("Precice configuration does not have 'FluidMesh'");
    }

    if (!_preciceInterface->hasMesh("BodyMesh")) {
      throwException("Precice configuration does not have 'BodyMesh'");
    }

    _fluidMeshId[0]
      = _preciceInterface->getMeshID("FluidMeshX");
    _fluidMeshId[1]
      = _preciceInterface->getMeshID("FluidMeshY");
    _bodyMeshId
      = _preciceInterface->getMeshID("BodyMesh");

    if (!_preciceInterface->hasData("VelocitiesX", _fluidMeshId[0])) {
      throwException("Precice configuration does not have 'Velocities' data"
                     " related to 'FluidMesh'");
    }

    if (!_preciceInterface->hasData("ForcesX", _fluidMeshId[0])) {
      throwException("Precice configuration does not have 'Forces' data"
                     " related to 'FluidMesh'");
    }

    _fluidMeshVelocitiesId[0]
      = _preciceInterface->getDataID("VelocitiesX", _fluidMeshId[0]);
    _fluidMeshVelocitiesId[1]
      = _preciceInterface->getDataID("VelocitiesY", _fluidMeshId[1]);
    _fluidMeshForcesId[0]
      = _preciceInterface->getDataID("ForcesX", _fluidMeshId[0]);
    _fluidMeshForcesId[1]
      = _preciceInterface->getDataID("ForcesY", _fluidMeshId[1]);
  }

  void
  computePositionInRespectToGeometry(CellAccessorType const& accessor) const {
    for (unsigned d = 0; d < Dimensions; ++d) {
      accessor.positionInRespectToGeometry()(d)
        = _preciceInterface->inquirePosition(
        accessor.velocityPosition(d).data(),
        std::set<int>({ _bodyMeshId }));
    }

    // if (!is_outside(accessor.positionInRespectToGeometry())) {
    // logInfo("{1}", accessor.pressurePosition().transpose());
    // }
  }

  void
  createFluidMeshVertex(CellAccessorType const& accessor) {
    for (unsigned d = 0; d < Dimensions; ++d) {
      unsigned distance
        = compute_cell_layer_along_geometry_interface(
        accessor,
        _maxLayerSize,
        d);

      bool doAdd = validate_layer_number(accessor,
                                         distance,
                                         outerLayerSize(),
                                         innerLayerSize(),
                                         d);

      if (doAdd) {
        // logInfo("{1} | d = {2}", accessor.index().transpose(), distance);
        VectorDsType position = accessor.velocityPosition(d);
        auto         vertexId
          = _preciceInterface->setMeshVertex(
          _fluidMeshId[d], position.data());

        this->_vertexIds[d].insert(
          std::make_pair(accessor.globalIndex(), vertexId));
      }
    }
  }

  std::pair<typename VertexIdMap::const_iterator, bool>
  doesVertexExist(CellAccessorType const& accessor,
                  unsigned const&         dimension) const {
    auto find_it = this->_vertexIds[dimension].find(accessor.globalIndex());

    if (find_it == this->_vertexIds[dimension].end()) {
      return std::make_pair(find_it, false);
    } else {
      return std::make_pair(find_it, true);
    }
  }

  void
  writeFluidVelocity(typename VertexIdMap::const_iterator const& it,
                     ScalarType const&                           velocity,
                     unsigned const&                             dimension) {
    _preciceInterface->writeScalarData(
      _fluidMeshVelocitiesId[dimension],
      it->second,
      static_cast<double>(velocity));
  }

  void
  mapData(ScalarType const& time_step_size) {
    _preciceInterface->advance(time_step_size);
  }

  void
  readFluidForce(typename VertexIdMap::const_iterator const& it,
                 ScalarType&                                 force,
                 unsigned const&                             dimension) {
    double temp;
    _preciceInterface->readScalarData(
      _fluidMeshForcesId[dimension],
      it->second,
      temp);
    force = static_cast<ScalarType>(temp);
  }

  VectorDsType
  computeBodyForceAt(unsigned          global_index,
                     ScalarType const& force,
                     MemoryType*       memory,
                     unsigned const&   dimension) {
    VectorDsType body_force = VectorDsType::Zero();

    auto accessor = *memory->grid()->begin();
    accessor.initialize(global_index);

    body_force(dimension)
      += accessor.width().prod() * (-force) / memory->timeStepSize();

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

  std::array<int, Dimensions> _fluidMeshId;
  int                         _bodyMeshId;
  std::array<int, Dimensions> _fluidMeshVelocitiesId;
  std::array<int, Dimensions> _fluidMeshForcesId;
};
}
}
}
