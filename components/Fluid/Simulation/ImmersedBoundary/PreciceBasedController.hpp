#pragma once

#include "Controller.hpp"

#include "Simulation/Configuration.hpp"
#include "Simulation/Grid.hpp"
#include "Simulation/IfsfdCellAccessor.hpp"
#include "Simulation/IfsfdMemory.hpp"
#include "Simulation/Private/mpigenerics.hpp"
#include "Simulation/SfsfdCellAccessor.hpp"
#include "Simulation/SfsfdMemory.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/ExecutionControl/exception>
#include <Uni/Helpers/macros>

#include <mpi.h>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

#include <map>
#include <vector>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
//

template <typename TVector>
class PreciceBasedInterfaceCell : public InterfaceCell<TVector> {
public:
  using Base = InterfaceCell<TVector>;

  enum {Dimensions = Base::Dimensions};

  using Vector = typename Base::Vector;

  using Scalar = typename Base::Scalar;

  using InternalIds = Eigen::Matrix<int, Dimensions, 1>;

  PreciceBasedInterfaceCell(unsigned const& global_index)
    : _globalIndex(global_index), _state(0) {}

  PreciceBasedInterfaceCell(unsigned const& global_index, int const& state)
    : _globalIndex(global_index), _state(state) {}

  unsigned
  globalIndex() const { return _globalIndex; }

  Vector const&
  data() const { return _data; }

  Vector&
  data() { return _data; }

  Scalar const&
  data(unsigned const& dimension) const {
    return _data(dimension);
  }

  Scalar&
  data(unsigned const& dimension) { return _data(dimension); }

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
using PreciceBasedInterfaceCellSet = boost::
                                     multi_index_container < std::shared_ptr
                                     < PreciceBasedInterfaceCell < TVector >>,
      boost::multi_index::
      indexed_by < boost::multi_index::
      ordered_unique < boost::multi_index::
      const_mem_fun < PreciceBasedInterfaceCell<TVector>,
      unsigned,
      &PreciceBasedInterfaceCell<TVector>::
      globalIndex >> >>;

template <typename TVector, typename TPredicate>
class PreciceBasedIteratorBackEnd :
  public IteratorBackEnd < InterfaceCell < TVector >> {
public:
  using Value = InterfaceCell<TVector>;

  using Base = IteratorBackEnd<Value>;

  using Set = PreciceBasedInterfaceCellSet<TVector>;

  using SetIterator = typename Set::template nth_index<0>::type::iterator;

  PreciceBasedIteratorBackEnd(SetIterator const& it,
                              SetIterator const& begin_it,
                              SetIterator const& end_it)
    : _it(it), _beginIt(begin_it), _endIt(end_it) {}

  PreciceBasedIteratorBackEnd(PreciceBasedIteratorBackEnd const&) = delete;

  ~PreciceBasedIteratorBackEnd() {}

  PreciceBasedIteratorBackEnd&
  operator=(PreciceBasedIteratorBackEnd const&)
    = delete;

  std::unique_ptr<Base> clone() const {
    return std::unique_ptr<Base>(
      new PreciceBasedIteratorBackEnd(_it, _beginIt, _endIt));
  }

  bool
  equals(std::unique_ptr<Base> const& other) {
    return _it
           == reinterpret_cast<PreciceBasedIteratorBackEnd const*>(other.get())
           ->_it;
  }

  Value&
  dereference() const { return *_it->get(); }

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

  using IteratorBackEndType = PreciceBasedIteratorBackEnd<Vector, TPredicate>;

  using Set = PreciceBasedInterfaceCellSet<TVector>;

  PreciceBasedIterableBackEnd(Set* set) : _set(set) {}

  PreciceBasedIterableBackEnd(PreciceBasedIterableBackEnd const&) = delete;

  ~PreciceBasedIterableBackEnd() {}

  PreciceBasedIterableBackEnd&
  operator=(PreciceBasedIterableBackEnd const&)
    = delete;

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
        return BaseIteratorType(BaseUniqueIteratorBackEndType(
                                  new IteratorBackEndType(it, it, it_end)));
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

    return BaseIteratorType(BaseUniqueIteratorBackEndType(
                              new IteratorBackEndType(it_end, it, it_end)));
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

    return BaseIteratorType(BaseUniqueIteratorBackEndType(
                              new IteratorBackEndType(find_it, it,
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

  enum {Dimensions = SolverTraitsType::Dimensions};

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

  using BaseIterableBackEndType = typename IterableType::BackEndType;

  using BaseSharedIterableBackEndType
          = typename IterableType::SharedBackEndType;

  struct VelocityPredicate;

  using VelocityIterableBackEndType
          = PreciceBasedIterableBackEnd<VectorDsType,
                                        typename PreciceBasedController::
                                        VelocityPredicate>;

  struct ForcePredicate;

  using ForceIterableBackEndType
          = PreciceBasedIterableBackEnd<VectorDsType,
                                        typename PreciceBasedController::
                                        ForcePredicate>;

  using Set = PreciceBasedInterfaceCellSet<VectorDsType>;

public:
  PreciceBasedController(Configuration const* configuration,
                         MemoryType const*    memory)
    : _memory(memory) {
    _outerLayerSize = configuration->get<unsigned>(
      "/Ib/Schemes/DirectForcing/PreciceBased/OuterLayerSize");
    _innerLayerSize = configuration->get<unsigned>(
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
    _fluidMeshForcesId = _preciceInterface->getDataID("Forces", _fluidMeshId);
  }

  void
  precompute() {
    // logInfo("Precice-based IB controller's precomutaions has started ...");
    _set.clear();
    _foreignCells.clear();
    // _preciceInterface->resetMesh(_fluidMeshId);

    for (auto const& accessor : _memory->grid()->innerGrid) {
      bool doAdd
        = validate_layer_number(accessor, _outerLayerSize, _innerLayerSize);

      if (doAdd) {
        // logInfo("{1} | d = {2}", accessor.index().transpose(), distance);

        _insertAndCreateVertex(
          accessor.globalIndex(), accessor.pressurePosition(), 3);
        _insertNeighbors(accessor);
      }
    }

    _locateForeignCells();
    _locateForeignCells();
  }

  void
  _insertAndCreateVertex(unsigned const&     global_index,
                         VectorDsType const& position,
                         int const&          state) {
    auto status = _set.emplace(new InterfaceCellType(global_index, state));

    if (status.second) {
      if (state > 1 || state < -1) {
        auto vertex_id = _preciceInterface->setMeshVertex(
          _fluidMeshId, position.template cast<double>().data());
        // logInfo("@@@ {1} | {2}",
        // _memory->parallelDistribution()->rank(),
        // position.transpose());
        (*status.first)->vertexId() = vertex_id;
      }
    } else {
      if ((*status.first)->state() == 1) {
        if (state > 1 || state == -2) {
          auto vertex_id = _preciceInterface->setMeshVertex(
            _fluidMeshId, position.template cast<double>().data());
          (*status.first)->vertexId() = vertex_id;
          (*status.first)->state(state);
          // logInfo("@!@ {1} | {2}",
          // _memory->parallelDistribution()->rank(),
          // position.transpose());
        }
      } else if ((*status.first)->state() == 2) {
        if (state > 2) {
          (*status.first)->state(state);
        }
      }
    }
  }

  void
  _insertNeighbors(CellAccessorType const& accessor) {
    for (unsigned d = 0; d < Dimensions; ++d) {
      auto     spatial_index = accessor.relativeIndex(d, -1);
      unsigned serial_index;

      if (_doesBelogsToThisSubdomain(spatial_index, serial_index, 1)) {
        _set.emplace(new InterfaceCellType(serial_index, 1));
      }

      spatial_index = accessor.relativeIndex(d, +1);
      auto position = accessor.pressurePosition(d, +1);

      if (_doesBelogsToThisSubdomain(spatial_index, serial_index, 2)) {
        _insertAndCreateVertex(serial_index, position, 2);
      }

      for (unsigned d2 = 0; d2 < Dimensions; ++d2) {
        spatial_index = accessor.relativeIndex(d, +1, d2, -1);

        if (_doesBelogsToThisSubdomain(spatial_index, serial_index, 1)) {
          _set.emplace(new InterfaceCellType(serial_index, 1));
        }
      }
    }
  }

  bool
  _doesBelogsToThisSubdomain(VectorDiType const& spatial_index,
                             unsigned&           serial_index,
                             int const&          state) {
    int rank;

    if (_memory->parallelDistribution()
        ->convertIndentedLocalIndexFromSpatialToSerial(
          spatial_index,
          rank,
          serial_index)) {
      return true;
    } else {
      auto status = _foreignCells.emplace(rank, Set());
      status.first->second.emplace(
        new InterfaceCellType(serial_index, state));

      return false;
    }
  }

  void
  _locateForeignCells() {
    struct GridTraits {
      using CellAccessorType
              = Uni::StructuredGrid::Basic::GlobalMultiIndex
                <void, void, Dimensions, GridTraits>;

      using GridType
              = Uni::StructuredGrid::Basic::Grid
                <void, void, Dimensions, GridTraits>;
    };

    typename GridTraits::GridType neighbors;
    neighbors.initialize(VectorDiType::Constant(3));

    auto corner
      = _memory->parallelDistribution()->index
        - VectorDiType::Constant(1);

    std::vector<unsigned> send_sizes;

    send_sizes.resize(neighbors.size().prod(), 0);
    receive_sizes.clear();
    receive_sizes.resize(neighbors.size().prod(), 0);

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> receive_requests;
    std::vector<MPI_Status>  statuses;

    send_requests.resize(2 * neighbors.size().prod());
    receive_requests.resize(2 * neighbors.size().prod());
    statuses.resize(4 * neighbors.size().prod());

    unsigned neighbors_size = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      if (rank < 0
          || rank == _memory->parallelDistribution()->rank()) {
        continue;
      }

      auto find_it = _foreignCells.find(rank);

      if (find_it != _foreignCells.end()) {
        send_sizes[neighbor.globalIndex()]
          = static_cast<unsigned>(find_it->second.size());
      }
      // logInfo("Send {1}", send_sizes[neighbor.globalIndex()]);
      MPI_Isend(&send_sizes[neighbor.globalIndex()],
                1,
                MPI_UNSIGNED,
                rank,
                1,
                _memory->parallelDistribution()->mpiCommunicator,
                &send_requests[neighbors_size]);
      MPI_Irecv(&receive_sizes[neighbor.globalIndex()],
                1,
                MPI_UNSIGNED,
                rank,
                1,
                _memory->parallelDistribution()->mpiCommunicator,
                &receive_requests[neighbors_size]);
      ++neighbors_size;
    }
    MPI_Waitall(neighbors_size,
                send_requests.data(),
                statuses.data());
    MPI_Waitall(neighbors_size,
                receive_requests.data(),
                statuses.data() +  neighbors.size().prod());

    std::vector < std::vector < unsigned >> send_serial_indices;
    std::vector < std::vector < int >> send_states;
    std::vector < std::vector < int >> receive_states;
    send_serial_indices.resize(_foreignCells.size());
    send_states.resize(_foreignCells.size());
    receive_serial_indices.clear();
    receive_serial_indices.resize(neighbors.size().prod());
    receive_states.resize(neighbors.size().prod());

    unsigned send_index    = 0;
    unsigned receive_index = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      auto find_it = _foreignCells.find(rank);

      if (find_it != _foreignCells.end()) {
        send_serial_indices[send_index].resize(find_it->second.size());
        send_states[send_index].resize(find_it->second.size());

        unsigned cell_index = 0;

        for (auto const& cell : find_it->second) {
          send_serial_indices[send_index][cell_index] = cell->globalIndex();
          send_states[send_index][cell_index]         = cell->state();
          ++cell_index;
        }

        MPI_Isend(
          send_serial_indices[send_index].data(),
          find_it->second.size(),
          MPI_UNSIGNED,
          rank,
          2,
          _memory->parallelDistribution()->mpiCommunicator,
          &send_requests[send_index]);
        MPI_Isend(
          send_states[send_index].data(),
          find_it->second.size(),
          MPI_INT,
          rank,
          3,
          _memory->parallelDistribution()->mpiCommunicator,
          &send_requests[send_index + neighbors.size().prod()]);
        ++send_index;
      }

      auto const temp_receive_size = receive_sizes[neighbor.globalIndex()];

      if (temp_receive_size == 0) {
        continue;
      }

      receive_states[neighbor.globalIndex()].resize(temp_receive_size);
      receive_serial_indices[neighbor.globalIndex()].resize(temp_receive_size);

      MPI_Irecv(
        receive_serial_indices[neighbor.globalIndex()].data(),
        temp_receive_size,
        MPI_UNSIGNED,
        rank,
        2,
        _memory->parallelDistribution()->mpiCommunicator,
        &receive_requests[receive_index]);
      MPI_Irecv(
        receive_states[neighbor.globalIndex()].data(),
        temp_receive_size,
        MPI_INT,
        rank,
        3,
        _memory->parallelDistribution()->mpiCommunicator,
        &receive_requests[receive_index + neighbors.size().prod()]);
      ++receive_index;
    }

    MPI_Waitall(send_index,
                send_requests.data(),
                statuses.data());
    MPI_Waitall(send_index,
                &send_requests[neighbors.size().prod()],
                &statuses[neighbors.size().prod()]);
    MPI_Waitall(receive_index,
                receive_requests.data(),
                &statuses[2 * neighbors.size().prod()]);
    MPI_Waitall(receive_index,
                &receive_requests[neighbors.size().prod()],
                &statuses[3 * neighbors.size().prod()]);

    send_size    = receive_index;
    receive_size = send_index;

    for (auto const& neighbor : neighbors) {
      auto const temp_receive_size = receive_sizes[neighbor.globalIndex()];

      for (unsigned i = 0; i < temp_receive_size; ++i) {
        auto serial_index = receive_serial_indices[neighbor.globalIndex()][i];

        auto accessor = *_memory->grid()->innerGrid.begin();
        accessor.initialize(serial_index);

        VectorDsType position = accessor.pressurePosition();

        _insertAndCreateVertex(serial_index,
                               position,
                               -receive_states[neighbor.globalIndex()][i]);

        for (unsigned d = 0; d < Dimensions; ++d) {
          auto spatial_index = accessor.relativeIndex(d, -1);

          if (_doesBelogsToThisSubdomain(spatial_index, serial_index, 1)) {
            _set.emplace(new InterfaceCellType(serial_index, 1));
          }
        }
      }
    }
  }

  void
  _exchangeData() {
    struct GridTraits {
      using CellAccessorType
              = Uni::StructuredGrid::Basic::GlobalMultiIndex
                <void, void, Dimensions, GridTraits>;

      using GridType
              = Uni::StructuredGrid::Basic::Grid
                <void, void, Dimensions, GridTraits>;
    };

    typename GridTraits::GridType neighbors;
    neighbors.initialize(VectorDiType::Constant(3));

    auto corner
      = _memory->parallelDistribution()->index - VectorDiType::Constant(1);

    std::vector < std::vector < ScalarType >> send_velocities;
    send_velocities.resize(send_size);

    std::vector < std::vector < ScalarType >> receive_velocities;
    receive_velocities.resize(receive_size);

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Status>  statuses;
    std::vector<MPI_Request> receive_requests;
    send_requests.resize(send_size);
    receive_requests.resize(receive_size);
    statuses.resize(send_size + receive_size);

    unsigned send_index    = 0;
    unsigned receive_index = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      auto find_it = _foreignCells.find(rank);

      if (find_it != _foreignCells.end()) {
        receive_velocities[receive_index].resize(Dimensions * find_it->second.size());

        MPI_Irecv(
          receive_velocities[receive_index].data(),
          Dimensions * find_it->second.size(),
          Private::getMpiScalarType<ScalarType>(),
          rank,
          4,
          _memory->parallelDistribution()->mpiCommunicator,
          &receive_requests[receive_index]);
        ++receive_index;
      }

      auto const temp_receive_size = receive_sizes[neighbor.globalIndex()];

      if (temp_receive_size == 0) {
        continue;
      }

      send_velocities[send_index].resize(Dimensions * temp_receive_size);

      for (unsigned i = 0; i < temp_receive_size; ++i) {
        auto find_it = _set.find(receive_serial_indices[neighbor.globalIndex()][i]);
        assert(find_it != _set.end());
        auto cell = *find_it;

        for (unsigned d = 0; d < Dimensions; ++d) {
          send_velocities[send_index][Dimensions * i + d] = cell->data(d);
        }
      }

      MPI_Isend(
        send_velocities[send_index].data(),
        Dimensions * temp_receive_size,
        Private::getMpiScalarType<ScalarType>(),
        rank,
        4,
        _memory->parallelDistribution()->mpiCommunicator,
        &send_requests[send_index]);
      ++send_index;
    }

    MPI_Waitall(send_size,
                send_requests.data(),
                statuses.data());
    MPI_Waitall(receive_size,
                receive_requests.data(),
                statuses.data() + send_size);

    receive_index = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      auto find_it = _foreignCells.find(rank);

      if (find_it != _foreignCells.end()) {
        unsigned cell_index = 0;

        for (auto const& cell : find_it->second) {
          for (unsigned d = 0; d < Dimensions; ++d) {
            cell->data(d) = receive_velocities[receive_index][cell_index + d];
          }
          // logInfo("from {1} to {2} | velocity = {3} | id = {4}",
          // _memory->parallelDistribution()->rank(),
          // rank,
          // cell->data(dimension),
          // cell->internalIds(dimension));
          cell_index += Dimensions;
        }
        ++receive_index;
      }
    }
  }

  void
  _exchangeForces() {
    struct GridTraits {
      using CellAccessorType
              = Uni::StructuredGrid::Basic::GlobalMultiIndex
                <void, void, Dimensions, GridTraits>;

      using GridType
              = Uni::StructuredGrid::Basic::Grid
                <void, void, Dimensions, GridTraits>;
    };

    typename GridTraits::GridType neighbors;
    neighbors.initialize(VectorDiType::Constant(3));

    auto corner
      = _memory->parallelDistribution()->index - VectorDiType::Constant(1);

    std::vector < std::vector < ScalarType >> send_velocities;
    send_velocities.resize(send_size);

    std::vector < std::vector < ScalarType >> receive_velocities;
    receive_velocities.resize(receive_size);

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Status>  statuses;
    std::vector<MPI_Request> receive_requests;
    send_requests.resize(send_size);
    receive_requests.resize(receive_size);
    statuses.resize(send_size + receive_size);

    unsigned send_index    = 0;
    unsigned receive_index = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      auto find_it = _foreignCells.find(rank);

      if (find_it != _foreignCells.end()) {
        receive_velocities[receive_index].resize(Dimensions * find_it->second.size());

        MPI_Irecv(
          receive_velocities[receive_index].data(),
          Dimensions * find_it->second.size(),
          Private::getMpiScalarType<ScalarType>(),
          rank,
          4,
          _memory->parallelDistribution()->mpiCommunicator,
          &receive_requests[receive_index]);
        ++receive_index;
      }

      auto const temp_receive_size = receive_sizes[neighbor.globalIndex()];

      if (temp_receive_size == 0) {
        continue;
      }

      send_velocities[send_index].resize(Dimensions * temp_receive_size);

      for (unsigned i = 0; i < temp_receive_size; ++i) {
        auto find_it = _set.find(receive_serial_indices[neighbor.globalIndex()][i]);
        assert(find_it != _set.end());
        auto cell = *find_it;

        for (unsigned d = 0; d < Dimensions; ++d) {
          send_velocities[send_index][Dimensions * i + d] = cell->force(d);
        }
      }

      MPI_Isend(
        send_velocities[send_index].data(),
        Dimensions * temp_receive_size,
        Private::getMpiScalarType<ScalarType>(),
        rank,
        4,
        _memory->parallelDistribution()->mpiCommunicator,
        &send_requests[send_index]);
      ++send_index;
    }

    MPI_Waitall(send_size,
                send_requests.data(),
                statuses.data());
    MPI_Waitall(receive_size,
                receive_requests.data(),
                statuses.data() + send_size);

    receive_index = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      auto find_it = _foreignCells.find(rank);

      if (find_it != _foreignCells.end()) {
        unsigned cell_index = 0;

        for (auto const& cell : find_it->second) {
          VectorDsType temp(&receive_velocities[receive_index][cell_index]);
          cell->force() = temp;
          cell_index   += Dimensions;
        }
        ++receive_index;
      }
    }
  }

  InterfaceCellType*
  _findCell(VectorDiType const& spatial_index) {
    unsigned serial_index;
    int      rank;

    if (_memory->parallelDistribution()
        ->convertIndentedLocalIndexFromSpatialToSerial(
          spatial_index,
          rank,
          serial_index)) {
      auto find_it = _set.template get<0>().find(serial_index);

      assert(find_it != _set.template get<0>().end());

      return find_it->get();
    } else {
      auto find_it = _foreignCells.find(rank);
      assert(find_it != _foreignCells.end());
      auto cell_it =   find_it->second.find(serial_index);
      assert(cell_it != find_it->second.end());

      return cell_it->get();
    }
  }

  void
  processVelocities() {
    _exchangeData();

    for (auto const& interface_cell : _set) {
      if (interface_cell->state() == 1
          || interface_cell->state() == -1) {
        continue;
      }

      VectorDsType resulting_velocity = interface_cell->data();

      for (unsigned d = 0; d < Dimensions; ++d) {
        auto accessor = *_memory->grid()->innerGrid.begin();
        accessor.initialize(interface_cell->globalIndex());
        auto spatial_neighbor_index = accessor.relativeIndex(d, -1);
        auto cell                   = _findCell(spatial_neighbor_index);
        resulting_velocity(d) += cell->data(d);
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
      if (interface_cell->state() == 1
          || interface_cell->state() == -1) {
        continue;
      }
      Eigen::Matrix<double, Dimensions, 1> temp;

      _preciceInterface->readVectorData(
        _fluidMeshForcesId, interface_cell->vertexId(), temp.data());
      interface_cell->force() = temp.template cast<ScalarType>();
    }

    _exchangeForces();

    for (auto& interface_cell : _set) {
      if (interface_cell->state() != 3) {
        continue;
      }

      VectorDsType force = interface_cell->force();

      for (unsigned d = 0; d < Dimensions; ++d) {
        auto accessor = *_memory->grid()->innerGrid.begin();
        accessor.initialize(interface_cell->globalIndex());

        auto spatial_neighbor_index = accessor.relativeIndex(d, +1);

        auto cell = _findCell(spatial_neighbor_index);

        force(d) += cell->force(d);
      }
      interface_cell->data() = force / (2.0 * _memory->timeStepSize());
    }
  }

  IterableType
  getVelocityIterable() {
    return IterableType(
      BaseSharedIterableBackEndType(new VelocityIterableBackEndType(&_set)));
  }

  IterableType
  getForceIterable() {
    return IterableType(
      BaseSharedIterableBackEndType(new ForceIterableBackEndType(&_set)));
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

  std::map<int, Set>    _foreignCells;
  std::vector<unsigned> receive_sizes;
  std::vector < std::vector < unsigned >> receive_serial_indices;
  unsigned send_size;
  unsigned receive_size;

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
