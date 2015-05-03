#pragma once

#include "Controller.hpp"

#include "Simulation/Configuration.hpp"
#include "Simulation/Grid.hpp"
#include "Simulation/IfsfdCellAccessor.hpp"
#include "Simulation/IfsfdMemory.hpp"
#include "Simulation/SfsfdCellAccessor.hpp"
#include "Simulation/SfsfdMemory.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/ExecutionControl/exception>
#include <Uni/Helpers/macros>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/multi_index/mem_fun.hpp>
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

  enum {
    Dimensions = Base::Dimensions
  };

  using Vector = typename Base::Vector;

  using Scalar = typename Base::Scalar;

  using InternalIds = Eigen::Matrix<int, Dimensions, 1>;

  PreciceBasedInterfaceCell(unsigned const& global_index)
    : _globalIndex(global_index) {}

  PreciceBasedInterfaceCell(unsigned const& global_index,
                            int const&      state)
    : _globalIndex(global_index),
    _state(state) {}

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

  Scalar const&
  force(unsigned const& dimension) const {
    return _force(dimension);
  }

private:
  unsigned _globalIndex;
  Uni_PrivateProperty(int,    vertexId);
  Uni_PrivateProperty(int,    state);
  Vector _data;
  Uni_PrivateProperty(Vector, force);
};

template <typename TVector>
using PreciceBasedInterfaceCellSet
        = boost::multi_index_container
          < std::shared_ptr < PreciceBasedInterfaceCell < TVector >>,
      boost::multi_index::indexed_by
      < boost::multi_index::ordered_unique
      < boost::multi_index::const_mem_fun
      < PreciceBasedInterfaceCell<TVector>,
      unsigned,
      &PreciceBasedInterfaceCell<TVector>::globalIndex >> >>;

template <typename TVector, typename TPredicate>
class PreciceBasedIteratorBackEnd :
  public IteratorBackEnd < InterfaceCell < TVector >> {
public:
  using Value = InterfaceCell<TVector>;

  using Base = IteratorBackEnd<Value>;

  using Set = PreciceBasedInterfaceCellSet<TVector>;

  using SetIterator
          = typename Set::template nth_index<0>::type::iterator;

  PreciceBasedIteratorBackEnd(SetIterator const& it,
                              SetIterator const& begin_it,
                              SetIterator const& end_it)
    : _it(it),
    _beginIt(begin_it),
    _endIt(end_it) {}

  PreciceBasedIteratorBackEnd(PreciceBasedIteratorBackEnd const&) = delete;

  ~PreciceBasedIteratorBackEnd() {}

  PreciceBasedIteratorBackEnd&
  operator=(PreciceBasedIteratorBackEnd const&) = delete;

  std::unique_ptr<Base>
  clone() const {
    return std::unique_ptr<Base>(
      new PreciceBasedIteratorBackEnd(_it,
                                      _beginIt,
                                      _endIt));
  }

  bool
  equals(std::unique_ptr<Base> const& other) {
    return _it
           == reinterpret_cast<PreciceBasedIteratorBackEnd const*>
           (other.get())->_it;
  }

  Value&
  dereference() const {
    return *_it->get();
  }

  void
  increment() {
    ++_it;

    for (; _it != _endIt; ++_it) {
      if (_predicate(*_it)) {
        return;
      }
    }
  }

  void
  decrement() {
    --_it;

    for (; _it != _beginIt; --_it) {
      if (_predicate(*_it)) {
        return;
      }
    }
  }

private:
  SetIterator _it;
  SetIterator _beginIt;
  SetIterator _endIt;
  TPredicate  _predicate;
};

template <typename TVector, typename TPredicate>
class PreciceBasedIterableBackEnd :
  public IterableBackEnd < InterfaceCell < TVector >> {
public:
  using Vector = TVector;

  using Value = InterfaceCell<Vector>;

  using Base = IterableBackEnd<Value>;

  using BaseIteratorType = typename Base::IteratorType;

  using BaseUniqueIteratorBackEndType
          = typename BaseIteratorType::UniqueBackEndType;

  using IteratorBackEndType
          = PreciceBasedIteratorBackEnd<Vector, TPredicate>;

  using Set = PreciceBasedInterfaceCellSet<TVector>;

  PreciceBasedIterableBackEnd(Set* set) : _set(set) {}

  PreciceBasedIterableBackEnd(PreciceBasedIterableBackEnd const&) = delete;

  ~PreciceBasedIterableBackEnd() {}

  PreciceBasedIterableBackEnd&
  operator=(PreciceBasedIterableBackEnd const&) = delete;

  unsigned
  size() {
    auto it     = _set->template get<0>().begin();
    auto it_end = _set->template get<0>().end();

    unsigned size = 0;

    for (; it != it_end; ++it) {
      if (_predicate(*it)) {
        ++size;
      }
    }

    return size;
  }

  BaseIteratorType
  begin() {
    auto it     = _set->template get<0>().begin();
    auto it_end = _set->template get<0>().end();

    for (; it != it_end; ++it) {
      if (_predicate(*it)) {
        return BaseIteratorType(
          BaseUniqueIteratorBackEndType(
            new IteratorBackEndType(
              it,
              it,
              it_end)));
      }
    }

    return end();
  }

  BaseIteratorType
  end() {
    auto it     = _set->template get<0>().begin();
    auto it_end = _set->template get<0>().end();

    for (; it != it_end; ++it) {
      if (_predicate(*it)) {
        break;
      }
    }

    return BaseIteratorType(
      BaseUniqueIteratorBackEndType(
        new IteratorBackEndType(
          it_end,
          it,
          it_end)));
  }

  BaseIteratorType
  find(unsigned const& global_index) {
    auto find_it = _set->template get<0>().find(global_index);

    if (find_it == _set->template get<0>().end()) {
      return end();
    }

    if (!_predicate(*find_it)) {
      return end();
    }

    auto it     = _set->template get<0>().begin();
    auto it_end = _set->template get<0>().end();

    for (; it != it_end; ++it) {
      if (_predicate(*it)) {
        break;
      }
    }

    return BaseIteratorType(
      BaseUniqueIteratorBackEndType(new IteratorBackEndType(
                                      find_it,
                                      it,
                                      _set->template get<0>().end())));
  }

private:
  TPredicate _predicate;
  Set*       _set;
};

template <typename TSolverTraits>
class PreciceBasedController :
  public Controller<typename TSolverTraits::VectorDsType> {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  using BaseType = Controller<VectorDsType>;

  using IterableType = typename BaseType::IterableType;

private:
  using BaseInterfaceCellType = typename BaseType::InterfaceCellType;

  using InterfaceCellType = PreciceBasedInterfaceCell<VectorDsType>;

  using BaseIterableBackEndType
          = typename IterableType::BackEndType;

  using BaseSharedIterableBackEndType
          = typename IterableType::SharedBackEndType;

  struct VelocityPredicate;

  using VelocityIterableBackEndType
          = PreciceBasedIterableBackEnd
            <VectorDsType, typename PreciceBasedController::VelocityPredicate>;

  struct ForcePredicate;

  using ForceIterableBackEndType
          = PreciceBasedIterableBackEnd
            <VectorDsType, typename PreciceBasedController::ForcePredicate>;

  using Set = PreciceBasedInterfaceCellSet<VectorDsType>;

public:
  PreciceBasedController(Configuration const* configuration,
                         MemoryType const*    memory) :
    _memory(memory) {
    _outerLayerSize
      = configuration->get<unsigned>(
      "/Ib/Schemes/DirectForcing/PreciceBased/OuterLayerSize");
    _innerLayerSize
      = configuration->get<unsigned>(
      "/Ib/Schemes/DirectForcing/PreciceBased/InnerLayerSize");
  }

  unsigned
  getMaxLayerSize() const {
    return std::max(_outerLayerSize, _innerLayerSize);
  }

  void
  initialize(precice::SolverInterface* precice_interface) {
    // logInfo("Controller is being initialized");
    //
    _preciceInterface = precice_interface;

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

        _insertAndCreateVertex(accessor.globalIndex(),
                               accessor.pressurePosition(),
                               3);
        _insertNeighbors(accessor);
      }
    }
  }

  void
  _insertAndCreateVertex(unsigned const&     global_index,
                         VectorDsType const& position,
                         int const&          state) {
    auto status = _set.emplace(
      new InterfaceCellType(global_index, state));

    if (status.second) {
      auto vertex_id = _preciceInterface->setMeshVertex(
        _fluidMeshId,
        position.template cast<double>().data());
      (*status.first)->vertexId() = vertex_id;
    } else {
      if ((*status.first)->state() == 1) {
        auto vertex_id = _preciceInterface->setMeshVertex(
          _fluidMeshId,
          position.template cast<double>().data());
        (*status.first)->vertexId() = vertex_id;
        (*status.first)->state(state);
      } else if ((*status.first)->state() == 2) {
        (*status.first)->state(state);
      }
    }
  }

  void
  _insertNeighbors(CellAccessorType const& accessor) {
    for (unsigned d = 0; d < Dimensions; ++d) {
      auto left_global_index = accessor.relativeGlobalIndex(d, -1);

      _set.emplace(
        new InterfaceCellType(left_global_index, 1));

      auto right_global_index = accessor.relativeGlobalIndex(d, +1);
      auto right_position     = accessor.pressurePosition(d, +1);
      _insertAndCreateVertex(right_global_index, right_position, 2);

      for (unsigned d2 = 0; d2 < Dimensions; ++d2) {
        auto left_global_index
          = accessor.relativeGlobalIndex(d, +1, d2, -1);

        _set.emplace(
          new InterfaceCellType(left_global_index, 1));
      }
    }
  }

  void
  processVelocities() {
    for (auto const& interface_cell : _set) {
      if (interface_cell->state() == 1) {
        continue;
      }

      VectorDsType resulting_velocity = interface_cell->data();

      for (unsigned d = 0; d < Dimensions; ++d) {
        auto accessor = *_memory->grid()->innerGrid.begin();
        accessor.initialize(interface_cell->globalIndex());
        auto neighbor_global_index = accessor.relativeGlobalIndex(d, -1);

        auto find_it
          = _set.template get<0>().find(neighbor_global_index);

        if (find_it == _set.template get<0>().end()) {
          throwException("Wow, this must never happen");
        } else {
          resulting_velocity(d) += (*find_it)->data(d);
        }
      }
      resulting_velocity /= 2.0;
      _preciceInterface->writeVectorData(
        _fluidMeshVelocitiesId,
        interface_cell->vertexId(),
        resulting_velocity.template cast<double>().data());
    }
  }

  void
  processForces() {
    for (auto& interface_cell : _set) {
      if (interface_cell->state() == 1) {
        continue;
      }
      Eigen::Matrix<double, Dimensions, 1> temp;

      _preciceInterface->readVectorData(
        _fluidMeshForcesId,
        interface_cell->vertexId(),
        temp.data());
      interface_cell->force() = temp.template cast<ScalarType>();
    }

    for (auto& interface_cell : _set) {
      if (interface_cell->state() != 3) {
        continue;
      }

      VectorDsType force = interface_cell->force();

      for (unsigned d = 0; d < Dimensions; ++d) {
        auto accessor = *_memory->grid()->innerGrid.begin();
        accessor.initialize(interface_cell->globalIndex());

        auto neighbor_global_index = accessor.relativeGlobalIndex(d, +1);

        auto find_it = _set.template get<0>().find(neighbor_global_index);

        if (find_it == _set.template get<0>().end()) {
          throwException("Wow, this must never happen");
        } else {
          force(d) += (*find_it)->force(d);
        }
      }
      interface_cell->data() = force / (2.0 * _memory->timeStepSize());
    }
  }

  IterableType
  getVelocityIterable() {
    return IterableType(
      BaseSharedIterableBackEndType(
        new VelocityIterableBackEndType(&_set)));
  }

  IterableType
  getForceIterable() {
    return IterableType(
      BaseSharedIterableBackEndType(
        new ForceIterableBackEndType(&_set)));
  }

private:
  precice::SolverInterface* _preciceInterface;
  MemoryType const*         _memory;
  unsigned                  _outerLayerSize;
  unsigned                  _innerLayerSize;

  int _fluidMeshId;
  int _bodyMeshId;
  int _fluidMeshVelocitiesId;
  int _fluidMeshForcesId;

  Set _set;

  struct VelocityPredicate {
    bool
    operator()(std::shared_ptr<InterfaceCellType> const& interface_cell) const {
      ((void)interface_cell);

      return true;
    }
  };

  struct ForcePredicate {
    bool
    operator()(std::shared_ptr<InterfaceCellType> const& interface_cell) const {
      if (interface_cell->state() == 3) {
        return true;
      }

      return false;
    }
  };
};
}
}
}
