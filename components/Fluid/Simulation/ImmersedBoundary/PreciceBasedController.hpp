#pragma once

#include "Controller.hpp"

#include "Simulation/Grid.hpp"
#include "Simulation/IfsfdCellAccessor.hpp"
#include "Simulation/IfsfdMemory.hpp"
#include "Simulation/SfsfdCellAccessor.hpp"
#include "Simulation/SfsfdMemory.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/ExecutionControl/exception>
#include <Uni/Helpers/macros>

// #include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
//

template <typename TVector>
class PreciceBasedInterfaceCell : public InterfaceCell<TVector> {
public:
  using Base = InterfaceCell<TVector>;

  using Vector = typename Base::Vector;

  using Scalar = typename Base::Scalar;

  PreciceBasedInterfaceCell(unsigned const& global_index,
                            int const&      vertexId)
    : _globalIndex(global_index),
    _vertexId(vertexId) {}

  unsigned
  globalIndex() const {
    return _globalIndex;
  }

  Vector const&
  data() const {
    return _data;
  }

  Vector&
  data() {
    return _data;
  }

  Scalar const&
  data(unsigned const& dimension) const {
    return _data(dimension);
  }

  Scalar&
  data(unsigned const& dimension) {
    return _data(dimension);
  }

  Uni_PublicProperty(int, vertexId);

private:
  unsigned _globalIndex;
  Vector   _data;
};

template <typename TVector>
using PreciceBasedInterfaceCellSet
        = boost::multi_index_container
          < PreciceBasedInterfaceCell<TVector>,
      boost::multi_index::indexed_by
      < boost::multi_index::ordered_unique
      < boost::multi_index::member
      < PreciceBasedInterfaceCell<TVector>,
      unsigned,
      &PreciceBasedInterfaceCell<TVector>::globalIndex >> >>;

template <typename TVector>
class PreciceBasedIteratorBackEnd :
  public IteratorBackEnd < InterfaceCell < TVector >> {
public:
  using Value = InterfaceCell<TVector>;

  using Base = IteratorBackEnd<Value>;

  using Set = PreciceBasedInterfaceCellSet<TVector>;

  using SetIterator = typename Set::template ntx_index<0>::type::iterator;

  PreciceBasedIteratorBackEnd(SetIterator it) : _it(it) {}

  virtual bool
  equals(Base const* other) {
    return _it == reinterpret_cast<PreciceBasedIteratorBackEnd*>(other)->_it;
  }

  virtual Value&
  dereference() const {
    return *_it;
  }

  virtual Value const&
  dereference() {
    return *_it;
  }

  virtual void
  increment() {
    ++_it;
  }

  virtual void
  decrement() {
    --_it;
  }

private:
  SetIterator _it;
};

template <typename TVector>
class PreciceBasedIterableBackEnd :
  public IterableBackEnd < InterfaceCell < TVector >> {
public:
  using Value = InterfaceCell<TVector>;

  using Base = IterableBackEnd<Value>;

  using IteratorType = typename Base::IteratorType;

  using Set = PreciceBasedInterfaceCellSet<TVector>;

  PreciceBasedIterableBackEnd(Set* set) : _set(set) {}

  unsigned
  size() {
    return static_cast<unsigned>(_set->size());
  }

  IteratorType
  begin() {
    return _set->begin();
  }

  IteratorType
  end() {
    return _set->end();
  }

  IteratorType
  find(unsigned const& global_index) {
    return _set->template get<0>().find(global_index);
  }

private:
  Set* _set;
};

template <typename TSolverTraits>
class PreciceBasedController :
  public Controller<typename TSolverTraits::VectorDsType> {
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

  using BaseType = Controller<VectorDsType>;

  using InterfaceCellType = typename BaseType::InterfaceCellType;

  using IterableType = typename BaseType::IterableType;

  using IterableBackEndType = PreciceBasedIterableBackEnd<VectorDsType>;

  using Set = PreciceBasedInterfaceCellSet<VectorDsType>;

public:
  PreciceBasedController(
    precice::SolverInterface* precice_interface,
    MemoryType const*         memory,
    unsigned const&           min_layer_size,
    unsigned const&           outer_layer_size,
    unsigned const&           inner_layer_size) :
    _preciceInterface(precice_interface),
    _memory(memory),
    _maxLayerSize(min_layer_size),
    _outerLayerSize(outer_layer_size),
    _innerLayerSize(inner_layer_size) {
    _maxLayerSize = std::max(_maxLayerSize, _outerLayerSize);
    _maxLayerSize = std::max(_maxLayerSize, _innerLayerSize);
  }

  void
  initialize() {
    // logInfo("Controller is being initialized");

    if (!_preciceInterface->hasMesh("FluidMesh")) {
      throwException("Precice configuration does not have 'FluidMesh'");
    }

    _fluidMeshId = _preciceInterface->getMeshID("FluidMesh");
    _bodyMeshId  = _preciceInterface->getMeshID("BodyMesh");

    if (!_preciceInterface->hasData("Velocities", _fluidMeshId)) {
      throwException("Precice configuration does not have 'Velocities' data"
                     " related to 'FluidMesh'");
    }

    if (!_preciceInterface->hasData("Forces", _fluidMeshId)) {
      throwException("Precice configuration does not have 'Forces' data"
                     " related to 'FluidMesh'");
    }

    _fluidMeshVelocitiesId
      = _preciceInterface->getDataID("Velocities", _fluidMeshId);
    _fluidMeshForcesId
      = _preciceInterface->getDataID("Forces", _fluidMeshId);
  }

  void
  precompute() {
    for (auto const& accessor : _memory->grid()->innerGrid) {
      bool doAdd = validate_layer_number(accessor,
                                         _outerLayerSize,
                                         _innerLayerSize);

      if (doAdd) {
        // logInfo("{1} | d = {2}", accessor.index().transpose(), distance);
        VectorDsType position = accessor.pressurePosition();

        auto vertex_id = _preciceInterface->setMeshVertex(
          _fluidMeshId, position.template cast<ScalarType>().data());

        _forceSet.emplace(accessor.globalindex(), vertex_id);
      }
    }
  }

  void
  processVelocities() {
    for (auto const& interface_cell : _forceSet) {
      VectorDsType resulting_velocity
        = interface_cell.data();

      for (unsigned d = 0; d < Dimensions; ++d) {
        auto accessor = *_memory->grid()->innerGrid.begin();
        accessor.initialize(interface_cell.globalIndex());
        auto neighbor_global_index = accessor.relativeGlobalIndex(d, -1);

        auto find_it
          = _forceSet.template get<0>().find(neighbor_global_index);

        if (find_it == _forceSet.template get<0>().end()) {
          // resulting_velocity(d) += 0.0;
        } else {
          resulting_velocity(d) += find_it->data(d);
        }
      }
      resulting_velocity /= 2.0;
      _preciceInterface->writeVectorData(
        _fluidMeshVelocitiesId,
        interface_cell.vertexId(),
        resulting_velocity.template cast<double>().data());
    }
  }

  void
  processForces() {
    this->resetIntermediateData();

    for (auto& interface_cell : _forceSet) {
      Eigen::Matrix<double, Dimensions, 1> temp;
      _preciceInterface->readVectorData(
        _fluidMeshForcesId,
        interface_cell.vertexId(),
        temp.data());
      interface_cell.data() = temp.template cast<ScalarType>();
    }

    for (auto& interface_cell : _forceSet) {
      VectorDsType force = interface_cell.data();

      for (unsigned d = 0; d < Dimensions; ++d) {
        auto accessor = *_memory->grid()->innerGrid.begin();
        accessor.initialize(interface_cell.globalIndex());

        auto neighbor_global_index = accessor.relativeGlobalIndex(d, +1);

        auto find_it = this->_vertexIds[0].find(neighbor_global_index);

        if (find_it == this->_vertexIds[0].end()) {
          // force(d) += 0.0;
        } else {
          force(d) += find_it->second->data(d);
        }
      }
      force /= 2.0;

      interface_cell.data() = force / _memory->timeStepSize();
    }
  }

  IterableType
  getVelocityIterable() {
    return IterableType(
      new CompoundIterableBackEnd<VectorDsType>(
        new IterableBackEndType(&_velocitySet),
        new IterableBackEndType(&_forceSet)));
  }

  IterableType
  getForceIterable() {
    return IterableType(new IterableBackEndType(&_forceSet));
  }

protected:
  precice::SolverInterface* _preciceInterface;
  MemoryType const*         _memory;
  unsigned                  _maxLayerSize;
  unsigned                  _outerLayerSize;
  unsigned                  _innerLayerSize;

  int _fluidMeshId;
  int _bodyMeshId;
  int _fluidMeshVelocitiesId;
  int _fluidMeshForcesId;

  Set _velocitySet;
  Set _forceSet;
};
}
}
}
