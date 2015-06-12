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

#include <Eigen/Dense>

#include <mpi.h>
#include <petscksp.h>

#include <Uni/ExecutionControl/exception>
#include <Uni/Helpers/macros>

#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

#include <algorithm>
#include <array>
#include <map>
#include <vector>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
//

template <typename TVector>
class RbfBasedInterfaceCell : public InterfaceCell<TVector> {
public:
  using Base = InterfaceCell<TVector>;

  enum {Dimensions = Base::Dimensions};

  using Vector = typename Base::Vector;

  using Scalar = typename Base::Scalar;

  using InternalIds = Eigen::Matrix<int, Dimensions, 1>;

  RbfBasedInterfaceCell(unsigned const& serial_index)
    : _globalIndex(serial_index),
    _internalIds(InternalIds::Constant(-1)) {}

  RbfBasedInterfaceCell(unsigned const& serial_index,
                        Vector const&   position)
    : _globalIndex(serial_index),
    _position(position),
    _internalIds(InternalIds::Constant(-1)) {}

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

  bool
  isInternal(unsigned const& dimension) const {
    return _internalIds(dimension) >= 0;
  }

  int const&
  internalIds(unsigned const& dimension) const {
    return _internalIds(dimension);
  }

  int&
  internalIds(unsigned const& dimension) {
    return _internalIds(dimension);
  }

private:
  unsigned _globalIndex;
  Uni_PrivateProperty(Vector,      position);
  Uni_PrivateProperty(InternalIds, internalIds);
  Vector _data;
};

template <typename TVector>
using RbfBasedInterfaceCellSet
        = boost::multi_index_container
          < std::shared_ptr < RbfBasedInterfaceCell < TVector >>,
      boost::multi_index::indexed_by
      < boost::multi_index::ordered_unique
      <
      boost::multi_index::const_mem_fun
      < RbfBasedInterfaceCell<TVector>,
      unsigned,
      &RbfBasedInterfaceCell<TVector>::globalIndex >> >>;

template <typename TVector>
class RbfBasedIteratorBackEnd : public IteratorBackEnd < InterfaceCell
                                < TVector >> {
public:
  using Value = InterfaceCell<TVector>;

  using Base = IteratorBackEnd<Value>;

  using Set = RbfBasedInterfaceCellSet<TVector>;

  using SetIterator = typename Set::template nth_index<0>::type::iterator;

  RbfBasedIteratorBackEnd(SetIterator const& it) : _it(it) {}

  RbfBasedIteratorBackEnd(RbfBasedIteratorBackEnd const&) = delete;

  ~RbfBasedIteratorBackEnd() {}

  RbfBasedIteratorBackEnd&
  operator=(RbfBasedIteratorBackEnd const&) = delete;

  std::unique_ptr<Base> clone() const {
    return std::unique_ptr<Base>(new RbfBasedIteratorBackEnd(_it));
  }

  bool
  equals(std::unique_ptr<Base> const& other) {
    return _it
           == reinterpret_cast<RbfBasedIteratorBackEnd const*>(other.get())
           ->_it;
  }

  Value&
  dereference() const { return *_it->get(); }

  void
  increment() { ++_it; }

  void
  decrement() { --_it; }

private:
  SetIterator _it;
};

template <typename TVector>
class RbfBasedIterableBackEnd : public IterableBackEnd < InterfaceCell
                                < TVector >> {
public:
  using Vector = TVector;

  using Value = InterfaceCell<Vector>;

  using Base = IterableBackEnd<Value>;

  using BaseIteratorType = typename Base::IteratorType;

  using BaseUniqueIteratorBackEndType
          = typename BaseIteratorType::UniqueBackEndType;

  using IteratorBackEndType = RbfBasedIteratorBackEnd<Vector>;

  using Set = RbfBasedInterfaceCellSet<TVector>;

  RbfBasedIterableBackEnd(Set* set) : _set(set) {}

  RbfBasedIterableBackEnd(RbfBasedIterableBackEnd const&) = delete;

  ~RbfBasedIterableBackEnd() {}

  RbfBasedIterableBackEnd&
  operator=(RbfBasedIterableBackEnd const&) = delete;

  unsigned
  size() { return static_cast<unsigned>(_set->size()); }

  BaseIteratorType
  begin() {
    return BaseIteratorType(BaseUniqueIteratorBackEndType(
                              new IteratorBackEndType(
                                _set->template get<0>().begin())));
  }

  BaseIteratorType
  end() {
    return BaseIteratorType(BaseUniqueIteratorBackEndType(
                              new IteratorBackEndType(
                                _set->template get<0>().end())));
  }

  BaseIteratorType
  find(unsigned const& serial_index) {
    return BaseIteratorType(BaseUniqueIteratorBackEndType(
                              new IteratorBackEndType(
                                _set->template get<0>().find(serial_index))));
  }

private:
  Set* _set;
};

template <typename TSolverTraits>
class RbfBasedController :
  public Controller<typename TSolverTraits::VectorDsType> {
public:
  class LagrangianNode;

public:
  using SolverTraitsType = TSolverTraits;

  enum {Dimensions = SolverTraitsType::Dimensions};

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

  using IterableType = typename BaseType::IterableType;

private:
  using BaseInterfaceCellType = typename BaseType::InterfaceCellType;

  using InterfaceCellType = RbfBasedInterfaceCell<VectorDsType>;

  using BaseIterableBackEndType = typename IterableType::BackEndType;

  using BaseSharedIterableBackEndType
          = typename IterableType::SharedBackEndType;

  using IterableBackEndType = RbfBasedIterableBackEnd<VectorDsType>;

  using Set = RbfBasedInterfaceCellSet<VectorDsType>;

public:
  RbfBasedController(Configuration const* configuration,
                     MemoryType const*    memory)
    : _memory(memory),
    _doComputeSupportRadius(true),
    _doComputeImqShape(true) {
    if (!configuration->isOfType<std::string>(
          "/Ib/Schemes/DirectForcing/RbfBased/SupportRadius")) {
      _doComputeSupportRadius = false;
      _suppportRadius         = static_cast<ScalarType>(configuration->get<long
                                                                           double>(
                                                          "/Ib/Schemes/DirectForcing/RbfBased/SupportRadius"));
    }

    if (!configuration->isOfType<std::string>(
          "/Ib/Schemes/DirectForcing/RbfBased/ImqShape")) {
      _doComputeImqShape = false;
      _imqShape          = static_cast<ScalarType>(configuration->get<long
                                                                      double>(
                                                     "/Ib/Schemes/DirectForcing/RbfBased/ImqShape"));
    }

    if (configuration->is("/Ib/Features/DevelopingStructure")) {
      _isDevelopingStructure = true;
    } else {
      _isDevelopingStructure = false;
    }
  }

  void
  initialize(precice::SolverInterface* precice_interface) {
    _preciceInterface = precice_interface;

    _bodyMeshSet.insert(_preciceInterface->getMeshID("BodyMesh"));

    if (_doComputeSupportRadius) {
      _suppportRadius = _memory->gridGeometry()->maxCellWidth().maxCoeff();
    }

    if (_doComputeImqShape) {
      _imqShape = _memory->parameters()->re()
                  / _memory->gridGeometry()->maxCellWidth().maxCoeff();
    }

    if (_isDevelopingStructure) {
      if (!_preciceInterface->hasData("Displacements", *_bodyMeshSet.begin())) {
        throwException("Precice configuration does not have 'Displacements' data"
                       " related to 'BodyMesh'");
      }

      _displacementsId
        = _preciceInterface->getDataID("Displacements", *_bodyMeshSet.begin());
    }

    logInfo("RBF shape = {1}",      _imqShape);
    logInfo("Support radius = {1}", _suppportRadius);
  }

  void
  precompute() {
    logInfo("Locate Interface Cells has been starting ...");

    // precice::MeshHandle const& mesh_handle =
    // _preciceInterface->getMeshHandle("BodyMesh");
    unsigned body_id = 0;

    auto     body_mesh_id     = _preciceInterface->getMeshID("BodyMesh");
    unsigned lagrangians_size = _preciceInterface->getMeshVertexSize(body_mesh_id);

    std::vector<int>    ids;
    std::vector<double> vertices;
    ids.resize(lagrangians_size);
    vertices.resize(Dimensions * lagrangians_size);

    for (unsigned i = 0; i < lagrangians_size; ++i) {
      ids[i] = i;
    }

    _preciceInterface->getMeshVertices(body_mesh_id,
                                       lagrangians_size,
                                       ids.data(),
                                       vertices.data());

    // double                               max_temp = -1;
    // double                               mix_temp = 100000000000;
    // Eigen::Matrix<double, Dimensions, 1> temp_coords2;

    _lagrangianNodes.clear();
    _lagrangianIds.clear();
    // _lagrangianNodes.reserve(mesh_handle.vertices().size());

    VectorDsType lower_boundary
      = _memory->gridGeometry()->minCellWidth().cwiseProduct(
      _memory->parallelDistribution()->corner.template cast<ScalarType>());
    VectorDsType upper_boundary
      =  lower_boundary
        + _memory->gridGeometry()->minCellWidth().cwiseProduct(
      _memory->parallelDistribution()->localCellSize.template cast<ScalarType>());

    // std::vector < Eigen::Matrix < double, Dimensions, 1 >> trash_coords;

    for (unsigned i = 0; i < ids.size(); ++i) {
      Eigen::Matrix<double, Dimensions, 1> temp_coords(&vertices[Dimensions * i]);

      // for (unsigned i = 0; i < trash_coords.size(); ++i) {
      // if (((trash_coords[i] - temp_coords).cwiseAbs().eval().array()
      // <= 10.0 * std::numeric_limits<double>::epsilon()).all()) {
      // logInfo("@@@@@@!!!!!!!!!!!!!!!!!!!!!!!!!! {1}",
      // trash_coords[i].transpose());
      // logInfo("@@@@@@!!!!!!!!!!!!!!!!!!!!!!!!!! {1}",
      // temp_coords.transpose());
      // logInfo("@@@@@@!!!!!!!!!!!!!!!!!!!!!!!!!! {1}", (trash_coords[i] -
      // temp_coords).eval().transpose());
      // throwException("");
      // }
      // }
      // trash_coords.emplace_back(temp_coords);

      // if (max_temp != -1) {
      // max_temp = std::max(max_temp, (temp_coords - temp_coords2).norm());
      // mix_temp = std::min(mix_temp, (temp_coords - temp_coords2).norm());
      // } else {
      // max_temp = 0;
      // }

      // temp_coords2 = temp_coords;

      VectorDsType lagrangian_cell_position
        = temp_coords.template cast<ScalarType>();

      if ((lagrangian_cell_position.array() <= upper_boundary.array()).all()
          && (lagrangian_cell_position.array() > lower_boundary.array()).all()) {
        _lagrangianNodes.emplace_back(body_id, lagrangian_cell_position);

        _lagrangianIds.emplace_back(ids[i]);

        ++body_id;
      }
    }

    if ((_memory->parallelDistribution()->rank() - 1) >= 0) {
      MPI_Status status;
      MPI_Recv(&_lagrangianIdOffset,
               1,
               MPI_UNSIGNED,
               (_memory->parallelDistribution()->rank() - 1),
               20,
               _memory->parallelDistribution()->mpiCommunicator,
               &status);
    } else {
      _lagrangianIdOffset = 0;
    }

    if ((_memory->parallelDistribution()->rank() + 1)
        < _memory->parallelDistribution()->rankSize()) {
      unsigned global_lagrangian_size = _lagrangianIdOffset + body_id;

      MPI_Send(&global_lagrangian_size,
               1,
               MPI_UNSIGNED,
               (_memory->parallelDistribution()->rank() + 1),
               20,
               _memory->parallelDistribution()->mpiCommunicator);
    }

    // logInfo("Lagrangian offset {1} {2} {3}",
    // _memory->parallelDistribution()->rank(),
    // _lagrangianIdOffset,
    // body_id);

    for (auto& node : _lagrangianNodes) {
      node.id() += _lagrangianIdOffset;
    }

    _velocitySet.clear();
    _forceSet.clear();

    for (unsigned d = 0; d < Dimensions; ++d) {
      _locateEulerianCells(d);
    }

    logInfo("Locate Interface Cells has finished");
  }

  void
  processVelocities() {}

  void
  processForces() {
    _lagrangianDisplacements.resize(Dimensions * _lagrangianIds.size(), 0.0);

    if (_isDevelopingStructure) {
      _preciceInterface->readBlockVectorData(_displacementsId,
                                             _lagrangianIds.size(),
                                             _lagrangianIds.data(),
                                             _lagrangianDisplacements.data());
    }

    for (unsigned d = 0; d < Dimensions; ++d) {
      _resolveForces(d);
    }
  }

  IterableType
  getVelocityIterable() {
    using CompoundIterableType = CompoundIterableBackEnd<BaseInterfaceCellType>;

    return IterableType(
      BaseSharedIterableBackEndType(
        new CompoundIterableType(
          typename CompoundIterableType::SharedIterableBackEndType(
            new IterableBackEndType(&_velocitySet)),
          typename CompoundIterableType::SharedIterableBackEndType(
            new IterableBackEndType(&_forceSet)))));
  }

  IterableType
  getForceIterable() {
    return IterableType(
      BaseSharedIterableBackEndType(new IterableBackEndType(&_forceSet)));
  }

private:
  void
  _locateEulerianCells(unsigned const& dimension) {
    _lagrangianNodesSupports[dimension].clear();
    _foreignEulerianCells[dimension].clear();

    auto         width = _memory->gridGeometry()->minCellWidth();
    VectorDsType grid_shift;
    grid_shift             = 0.5 * width;
    grid_shift(dimension) += 0.5 * width(dimension);

    VectorDiType support_size
      = (VectorDsType::Constant(_suppportRadius).cwiseQuotient(width)
         + VectorDsType::Ones()).template cast<int>();

    unsigned internal_id = 0;

    logInfo("Size {1} ", _lagrangianNodes.size());

    for (auto& lagrangian_node : _lagrangianNodes) {
      _lagrangianNodesSupports[dimension].emplace_back(&lagrangian_node);

      VectorDiType size = (lagrangian_node.position() - grid_shift)
                          .cwiseQuotient(width)
                          .template cast<int>();

      // Compute search boundaries for the fluid cells
      // around the body vertex

      VectorDiType search_offset = size - support_size
                                   - 0 * VectorDiType::Ones();
      VectorDiType search_size = size + support_size + 2 * VectorDiType::Ones();

      struct GridTraits {
        using CellAccessorType
                = Uni::StructuredGrid::Basic::GlobalMultiIndex
                  <void, void, Dimensions, GridTraits>;

        using GridType
                = Uni::StructuredGrid::Basic::Grid
                  <void, void, Dimensions, GridTraits>;
      };

      typename GridTraits::GridType grid;
      grid.initialize(search_size - search_offset);

      for (auto const& accessor : grid) {
        VectorDiType global_spatial_index = search_offset + accessor.index();

        VectorDsType position
          = global_spatial_index.template cast<ScalarType>().cwiseProduct(width)
            + grid_shift;

        auto dist = (lagrangian_node.position() - position).norm();

        if (dist > _suppportRadius) {
          continue;
        }

        // Test if the eulerian cell belongs to this domain
        int                rank;
        unsigned           serial_index;
        InterfaceCellType* eulerian_cell;

        if (_memory->parallelDistribution()
            ->convertUnindentedGlobalIndexFromSpatialToSerial(
              global_spatial_index,
              rank,
              serial_index)) {
          eulerian_cell = _insertEulerianCell(dimension,
                                              serial_index,
                                              internal_id,
                                              position);
        } else {
          eulerian_cell = _insertReceiveForeignEulerianCell(dimension,
                                                            rank,
                                                            serial_index,
                                                            position);
        }

        _lagrangianNodesSupports[dimension].back().addEulerianCell(eulerian_cell);
      }
    }

    _locateForeignCells(dimension, internal_id);

    _localEulerianCellSize(dimension) = internal_id;

    if ((_memory->parallelDistribution()->rank() - 1) >= 0) {
      MPI_Status status;
      MPI_Recv(&_eulerianIdOffset[dimension],
               1,
               MPI_UNSIGNED,
               (_memory->parallelDistribution()->rank() - 1),
               21 + dimension,
               _memory->parallelDistribution()->mpiCommunicator,
               &status);
    } else {
      _eulerianIdOffset[dimension] = 0;
    }

    if ((_memory->parallelDistribution()->rank() + 1)
        < _memory->parallelDistribution()->rankSize()) {
      unsigned global_lagrangian_size = _eulerianIdOffset[dimension] + internal_id;

      MPI_Send(&global_lagrangian_size,
               1,
               MPI_UNSIGNED,
               (_memory->parallelDistribution()->rank() + 1),
               21 + dimension,
               _memory->parallelDistribution()->mpiCommunicator);
    }

    for (auto& cell : _forceSet) {
      if (cell->internalIds(dimension) != -1) {
        cell->internalIds(dimension) += _eulerianIdOffset[dimension];
      }
    }

    for (auto& cell : _velocitySet) {
      if (cell->internalIds(dimension) != -1) {
        cell->internalIds(dimension) += _eulerianIdOffset[dimension];
      }
    }
  }

  void
  _locateForeignCells(unsigned const& dimension,
                      unsigned&       internal_id) {
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

    auto corner = _memory->parallelDistribution()->index - VectorDiType::Constant(1);

    std::vector<unsigned> send_sizes;

    send_sizes.resize(neighbors.size().prod(), 0);
    receive_sizes[dimension].clear();
    receive_sizes[dimension].resize(neighbors.size().prod(), 0);

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> receive_requests;
    std::vector<MPI_Status>  statuses;

    send_requests.resize(neighbors.size().prod());
    receive_requests.resize(neighbors.size().prod());
    statuses.resize(2 * neighbors.size().prod());

    unsigned neighbors_size = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      if (rank < 0
          || rank == _memory->parallelDistribution()->rank()) {
        continue;
      }

      auto find_it = _foreignEulerianCells[dimension].find(rank);

      if (find_it != _foreignEulerianCells[dimension].end()) {
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
      MPI_Irecv(&receive_sizes[dimension][neighbor.globalIndex()],
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
    send_serial_indices.resize(_foreignEulerianCells[dimension].size());
    receive_serial_indices[dimension].clear();
    receive_serial_indices[dimension].resize(neighbors.size().prod());

    unsigned send_index    = 0;
    unsigned receive_index = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      auto find_it = _foreignEulerianCells[dimension].find(rank);

      if (find_it != _foreignEulerianCells[dimension].end()) {
        send_serial_indices[send_index].resize(find_it->second.size());

        unsigned cell_index = 0;

        for (auto const& cell : find_it->second) {
          send_serial_indices[send_index][cell_index] = cell->globalIndex();
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
        ++send_index;
      }

      auto const temp_receive_size
        = receive_sizes[dimension][neighbor.globalIndex()];

      if (temp_receive_size == 0) {
        continue;
      }

      receive_serial_indices[dimension][neighbor.globalIndex()]
      .resize(temp_receive_size);

      MPI_Irecv(
        receive_serial_indices[dimension][neighbor.globalIndex()].data(),
        temp_receive_size,
        MPI_UNSIGNED,
        rank,
        2,
        _memory->parallelDistribution()->mpiCommunicator,
        &receive_requests[receive_index]);
      ++receive_index;
    }

    MPI_Waitall(send_index,
                send_requests.data(),
                statuses.data());
    MPI_Waitall(receive_index,
                receive_requests.data(),
                &statuses[neighbors.size().prod()]);

    send_size[dimension]    = receive_index;
    receive_size[dimension] = send_index;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      auto const temp_receive_size
        = receive_sizes[dimension][neighbor.globalIndex()];

      for (unsigned i = 0; i < temp_receive_size; ++i) {
        _insertSendForeignEulerianCell(
          dimension,
          rank,
          receive_serial_indices[dimension][neighbor.globalIndex()][i],
          internal_id);
      }
    }
  }

  InterfaceCellType*
  _insertEulerianCell(unsigned const&     dimension,
                      unsigned const&     serial_index,
                      unsigned&           internal_id,
                      VectorDsType const& position) {
    // if (static_cast<int>(serial_index)
    // < _memory->grid()->innerGrid.begin()->globalIndex()
    // || static_cast<int>(serial_index)
    // > (_memory->grid()->innerGrid.end())->globalIndex()) {
    // logInfo("{1}", serial_index);
    // logInfo("{1}", _memory->grid()->innerGrid.begin()->globalIndex());
    // logInfo("{1}", (_memory->grid()->innerGrid.end())->globalIndex());
    // }

    assert(static_cast<int>(serial_index)
           >= _memory->grid()->innerGrid.begin()->globalIndex()
           && static_cast<int>(serial_index)
           <= (--_memory->grid()->innerGrid.end())->globalIndex());

    bool isInternal
      = (_preciceInterface->inquirePosition(
           position.template cast<double>().data(), _bodyMeshSet)
         != precice::constants::positionOutsideOfGeometry());

    auto find_it = _velocitySet.template get<0>().find(serial_index);

    if (find_it == _velocitySet.template get<0>().end()) {
      find_it = _forceSet.template get<0>().find(serial_index);

      if (find_it == _forceSet.template get<0>().end()) {
        if (isInternal) {
          auto status = _forceSet.emplace(
            std::make_shared<InterfaceCellType>(serial_index, position));
          (*status.first)->internalIds(dimension) = internal_id;
          ++internal_id;

          return status.first->get();
        } else {
          auto status = _velocitySet.emplace(
            std::make_shared<InterfaceCellType>(serial_index, position));

          return status.first->get();
        }
      } else {
        if (isInternal) {
          if ((*find_it)->internalIds(dimension) == -1) {
            (*find_it)->internalIds(dimension) = internal_id;
            ++internal_id;
          }
        }

        return find_it->get();
      }
    } else {
      if (isInternal) {
        auto status = _forceSet.emplace(*find_it);
        _velocitySet.erase(find_it);
        (*status.first)->internalIds(dimension) = internal_id;
        ++internal_id;

        return status.first->get();
      } else {
        return find_it->get();
      }
    }
  }

  InterfaceCellType*
  _insertReceiveForeignEulerianCell(unsigned const&     dimension,
                                    int const&          rank,
                                    unsigned const&     serial_index,
                                    VectorDsType const& position) {
    auto status = _foreignEulerianCells[dimension].emplace(
      std::make_pair(rank, Set()));

    auto cell = std::make_shared<InterfaceCellType>(serial_index, position);

    auto cell_status = status.first->second.template get<0>().emplace(cell);

    return cell_status.first->get();
  }

  InterfaceCellType*
  _insertSendForeignEulerianCell(unsigned const& dimension,
                                 int const&      rank,
                                 unsigned const& serial_index,
                                 unsigned&       internal_id) {
    auto spatial_index
      = _memory->parallelDistribution()
        ->convertSerialIndexToUnindentedGlobal(rank,
                                               serial_index);

    VectorDsType position
      = _memory->gridGeometry()->minCellWidth().cwiseProduct(
      spatial_index.template cast<ScalarType>())
        + 0.5 * _memory->gridGeometry()->minCellWidth();
    position(dimension)
      += 0.5 * _memory->gridGeometry()->minCellWidth(dimension);

    return _insertEulerianCell(dimension,
                               serial_index,
                               internal_id,
                               position);
  }

  InterfaceCellType*
  _findSendForeignEulerianCell(unsigned const& serial_index) {
    auto find_it = _velocitySet.template get<0>().find(serial_index);

    if (find_it == _velocitySet.template get<0>().end()) {
      find_it = _forceSet.template get<0>().find(serial_index);

      if (find_it == _forceSet.template get<0>().end()) {
        throwException("Wow this must never happen");
      }
    }

    return find_it->get();
  }

  void
  _resolveForces(unsigned const& dimension) {
    using Matrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

    // unsigned max_size = 0;
    // unsigned mix_size = 1000000000;

    for (auto& supports : _lagrangianNodesSupports[dimension]) {
      Matrix m(supports.size() + 1, supports.size() + 1);

      // max_size = std::max(max_size, static_cast<unsigned>(supports.size()));
      // mix_size = std::min(mix_size, static_cast<unsigned>(supports.size()));

      Vector b(supports.size() + 1);

      for (unsigned j = 0; j < supports.size(); ++j) {
        auto const& fluid_cell = supports[j];

        for (unsigned k = 0; k < j; ++k) {
          auto const& fluid_cell2 = supports[k];
          m(j, k) = computeImqRbf(fluid_cell.cell()->position(),
                                  fluid_cell2.cell()->position());
        }

        for (unsigned k = j + 1; k < supports.size(); ++k) {
          auto const& fluid_cell2 = supports[k];
          m(j, k) = computeImqRbf(fluid_cell.cell()->position(),
                                  fluid_cell2.cell()->position());
        }

        m(j, j)               = computeImqRbf(0.0);
        m(j, supports.size()) = 1.0;
        b(j)                  = computeImqRbf(
          supports.lagrangianNode()->position(),
          fluid_cell.cell()->position());
      }

      m.row(supports.size()).fill(1.0);
      m(supports.size(), supports.size()) = 0.0;
      b(supports.size())                  = 1.0;

      Eigen::ColPivHouseholderQR<Matrix> dec(m);
      Vector                             x = dec.solve(b);

      for (unsigned j = 0; j < supports.size(); ++j) {
        supports[j].weight() = x(j);
      }
    }

    // logInfo("Max in {1} is {2}", dimension, max_size);
    // logInfo("Mix in {1} is {2}", dimension, mix_size);

    _performMpiExchange(dimension);

    /*
     * Solve least-square problem
     */

    // Create right-hand side vector and vector of unknowns
    // logInfo("Create b = {1}", _lagrangianNodes.size());
    Vec b, x;
    Vec bt;
    VecCreate(_memory->parallelDistribution()->mpiCommunicator, &b);
    VecSetSizes(b, _lagrangianNodes.size(), PETSC_DECIDE);
    VecSetType(b, VECMPI);
    VecCreate(_memory->parallelDistribution()->mpiCommunicator, &bt);
    VecSetSizes(bt, _localEulerianCellSize(dimension), PETSC_DECIDE);
    VecSetType(bt, VECMPI);
    // logInfo("Create on {1} x = {2}",
    // _memory->parallelDistribution()->rank(),
    // _localEulerianCellSize(dimension));
    VecCreate(_memory->parallelDistribution()->mpiCommunicator, &x);
    VecSetSizes(x, _localEulerianCellSize(dimension), PETSC_DECIDE);
    VecSetType(x, VECMPI);

    // Fill right-hand side;

    // double temp_row         = 0;
    // double temp_norm_base   = 0;
    // double temp_norm_base_b = 0;

    unsigned index = 0;

    for (auto& supports : _lagrangianNodesSupports[dimension]) {
      PetscScalar value
        = _lagrangianDisplacements[Dimensions * index + dimension]
          / _memory->timeStepSize();

      // logInfo("{1}", _lagrangianDisplacements[Dimensions * index +
      // dimension]);

      for (auto& fluid_cell : supports) {
        value -= fluid_cell.weight() * fluid_cell.cell()->data(dimension);

        // temp_norm_base += fluid_cell.cell()->data(dimension)
        /// _memory->timeStepSize();

        // temp_norm_base_b += fluid_cell.weight()
        // * fluid_cell.cell()->data(dimension)
        /// _memory->timeStepSize();

        // if (supports.lagrangianNode()->id() == 1) {
        // temp_row += fluid_cell.weight();
        // }
      }

      value /= _memory->timeStepSize();

      PetscInt serial_index = supports.lagrangianNode()->id();
      // logInfo("b({1}) = {2}", serial_index, value);
      VecSetValues(b, 1, &serial_index, &value, INSERT_VALUES);
      ++index;
    }

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    // logInfo("========= ROW {1}", temp_row);

    // Create matrix
    // logInfo("Create matrix");
    Mat global_matrix;
    Mat preconditioning_global_matrix;
    MatCreate(_memory->parallelDistribution()->mpiCommunicator, &global_matrix);
    MatSetSizes(global_matrix,
                _lagrangianNodes.size(),
                _localEulerianCellSize(dimension),
                PETSC_DECIDE,
                PETSC_DECIDE);

    // MatSetType(global_matrix, MATMPIAIJ);
    MatSetType(global_matrix, MATAIJ);
    MatSetUp(global_matrix);

    // Fill matrix
    // PetscInt end;
    // MatGetOwnershipRange(global_matrix, &first, &end);
    // logInfo("Fill matrix");

    // double temp_column = 0;

    for (auto const& supports : _lagrangianNodesSupports[dimension]) {
      PetscInt global_row = supports.lagrangianNode()->id();

      for (auto const& fluid_cell : supports) {
        PetscInt global_column = fluid_cell.cell()->internalIds(dimension);

        if (!fluid_cell.cell()->isInternal(dimension)) {
          continue;
        }

        PetscScalar value = fluid_cell.weight();

        // logInfo("m({1}, {2}) = {3}", global_row, global_column, value);

        MatSetValues(global_matrix,
                     1,
                     &global_row,
                     1,
                     &global_column,
                     &value,
                     INSERT_VALUES);

        // if (global_column == 1) {
        // temp_column += value;
        // }
      }
    }

    MatAssemblyBegin(global_matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(global_matrix, MAT_FINAL_ASSEMBLY);

    // logInfo("========= COLUMN {1}", temp_column);

    KSP solver;
    KSPCreate(_memory->parallelDistribution()->mpiCommunicator, &solver);
    KSPSetType(solver, KSPFGMRES);
    PC preconditioner;
    KSPGetPC(solver, &preconditioner);
    PCSetType(preconditioner, PCKSP);
    KSP preconditioner_solver;
    PCKSPGetKSP(preconditioner, &preconditioner_solver);
    KSPSetType(preconditioner_solver, KSPGMRES);
    KSPSetTolerances(preconditioner_solver, 1e-21, 1e-25, PETSC_DEFAULT, 10);
    // PCSetType(preconditioner, PCILU);
    KSPSetTolerances(solver,                1e-21, 1e-25, PETSC_DEFAULT, 10000);
    MatTransposeMatMult(global_matrix,
                        global_matrix,
                        MAT_INITIAL_MATRIX,
                        PETSC_DEFAULT,
                        &preconditioning_global_matrix);
#if ((PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 5))
    KSPSetOperators(solver,
                    // global_matrix,
                    preconditioning_global_matrix,
                    preconditioning_global_matrix
                    // 0
                    );
#else
    KSPSetOperators(solver,
                    // global_matrix,
                    preconditioning_global_matrix,
                    preconditioning_global_matrix,
                    // 0,
                    SAME_PRECONDITIONER);
#endif

    // logInfo("KSP Setup");
    KSPSetUp(solver);
    // logInfo("KSP Solve");
    MatMultTranspose(global_matrix, b, bt);
    KSPSolve(solver, bt, x);

    // MatMult(global_matrix, x, b);
    // PetscScalar temp_norm;
    // VecNorm(x, NORM_1, &temp_norm);
    // logInfo("NOOOOOOOOOOOOORM1 {1}", temp_norm + temp_norm_base);
    // VecNorm(b, NORM_1, &temp_norm);
    // logInfo("NOOOOOOOOOOOOORM2 {1}", temp_norm + temp_norm_base_b);

#if !defined (NDEBUG)
    PetscScalar global_norm;
    PetscInt    iteration_number;
    VecNorm(x, NORM_2, &global_norm);
    KSPGetIterationNumber(solver, &iteration_number);

    if (_memory->parallelDistribution()->rank() == 0) {
      logInfo("Norm of vector {1} iterations {2}", global_norm, iteration_number);
    }
#endif

    // Retrieve data
    // logInfo("Retrieve data");
    PetscScalar* unknowns;
    VecGetArray(x, &unknowns);

    for (auto& cell : _forceSet) {
      if (cell->isInternal(dimension)) {
        // if ((_memory->parallelDistribution()->rank()) == 0) {
        // logInfo("{1} {2}",
        // unknowns[cell->internalIds(dimension) -
        // _eulerianIdOffset[dimension]],
        // cell->internalIds(dimension));
        // }
        cell->data(dimension)
          = unknowns[cell->internalIds(dimension) - _eulerianIdOffset[dimension]];
      } else {
        cell->data(dimension) = 0.0;
      }
    }

    VecRestoreArray(x, &unknowns);

    // logInfo("Clean up");
    VecDestroy(&b);
    VecDestroy(&bt);
    VecDestroy(&x);
    MatDestroy(&global_matrix);
    MatDestroy(&preconditioning_global_matrix);
    KSPDestroy(&solver);
  }

  void
  _performMpiExchange(unsigned const& dimension) {
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

    std::vector < std::vector < ScalarType >> send_velocities;
    std::vector < std::vector < int >> send_internal_ids;

    send_velocities.resize(send_size[dimension]);
    send_internal_ids.resize(send_size[dimension]);

    std::vector < std::vector < ScalarType >> receive_velocities;
    std::vector < std::vector < int >> receive_internal_ids;

    receive_velocities.resize(receive_size[dimension]);
    receive_internal_ids.resize(receive_size[dimension]);

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Status>  statuses;
    std::vector<MPI_Request> receive_requests;
    send_requests.resize(2 * send_size[dimension]);
    receive_requests.resize(2 * receive_size[dimension]);
    statuses.resize(2 * send_size[dimension] + 2 * receive_size[dimension]);

    unsigned send_index    = 0;
    unsigned receive_index = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      auto find_it = _foreignEulerianCells[dimension].find(rank);

      if (find_it != _foreignEulerianCells[dimension].end()) {
        receive_velocities[receive_index].resize(find_it->second.size());
        receive_internal_ids[receive_index].resize(find_it->second.size());

        MPI_Irecv(
          receive_velocities[receive_index].data(),
          find_it->second.size(),
          Private::getMpiScalarType<ScalarType>(),
          rank,
          3,
          _memory->parallelDistribution()->mpiCommunicator,
          &receive_requests[receive_index]);

        MPI_Irecv(
          receive_internal_ids[receive_index].data(),
          find_it->second.size(),
          MPI_INT,
          rank,
          4,
          _memory->parallelDistribution()->mpiCommunicator,
          &receive_requests[receive_index +  receive_size[dimension]]);
        ++receive_index;
      }

      auto const temp_receive_size
        = receive_sizes[dimension][neighbor.globalIndex()];

      if (temp_receive_size == 0) {
        continue;
      }

      send_velocities[send_index].resize(temp_receive_size);
      send_internal_ids[send_index].resize(temp_receive_size);

      for (unsigned i = 0; i < temp_receive_size; ++i) {
        auto cell = _findSendForeignEulerianCell(
          receive_serial_indices[dimension][neighbor.globalIndex()][i]);

        send_velocities[send_index][i]   = cell->data(dimension);
        send_internal_ids[send_index][i] = cell->internalIds(dimension);
      }

      MPI_Isend(
        send_velocities[send_index].data(),
        temp_receive_size,
        Private::getMpiScalarType<ScalarType>(),
        rank,
        3,
        _memory->parallelDistribution()->mpiCommunicator,
        &send_requests[send_index]);

      MPI_Isend(
        send_internal_ids[send_index].data(),
        temp_receive_size,
        MPI_INT,
        rank,
        4,
        _memory->parallelDistribution()->mpiCommunicator,
        &send_requests[send_index + send_size[dimension]]);
      ++send_index;
    }

    MPI_Waitall(send_size[dimension],
                send_requests.data(),
                statuses.data());
    MPI_Waitall(send_size[dimension],
                send_requests.data() + send_size[dimension],
                statuses.data() + send_size[dimension]);
    MPI_Waitall(receive_size[dimension],
                receive_requests.data(),
                statuses.data() + 2 * send_size[dimension]);
    MPI_Waitall(receive_size[dimension],
                receive_requests.data() + receive_size[dimension],
                statuses.data()
                + 2 * send_size[dimension] + receive_size[dimension]);

    receive_index = 0;

    for (auto const& neighbor : neighbors) {
      int rank = _memory->parallelDistribution()
                 ->getMpiRankFromSubdomainSpatialIndex(
        corner + neighbor.index());

      auto find_it = _foreignEulerianCells[dimension].find(rank);

      if (find_it != _foreignEulerianCells[dimension].end()) {
        unsigned cell_index = 0;

        for (auto const& cell : find_it->second) {
          cell->data(dimension)
            = receive_velocities[receive_index][cell_index];
          cell->internalIds(dimension)
            = receive_internal_ids[receive_index][cell_index];
          // logInfo("from {1} to {2} | velocity = {3} | id = {4}",
          // _memory->parallelDistribution()->rank(),
          // rank,
          // cell->data(dimension),
          // cell->internalIds(dimension));
          ++cell_index;
        }
        ++receive_index;
      }
    }
  }

  //

public:
  class LagrangianNode {
public:
    LagrangianNode(unsigned const& id, VectorDsType const& position)
      : _id(id), _position(position) {}

    Uni_PublicProperty(unsigned,     id);
    Uni_PublicProperty(VectorDsType, position);
  };

  class LagrangianNodeSupports {
public:
    class Support {
public:
      Support(InterfaceCellType const* cell, ScalarType const& weight)
        : _cell(cell), _weight(weight) {}

      Uni_PublicProperty(InterfaceCellType const*, cell);
      Uni_PublicProperty(ScalarType,               weight);
    };
    using Supports = std::vector<Support>;

    LagrangianNodeSupports(LagrangianNode* lagrangian_node)
      : _lagrangianNode(lagrangian_node) {}

    typename Supports::const_iterator
    begin() const {
      return _supports.begin();
    }

    typename Supports::iterator
    begin() { return _supports.begin(); }

    typename Supports::const_iterator
    end() const { return _supports.end(); }

    typename Supports::iterator
    end() { return _supports.end(); }

    std::size_t
    size() const { return _supports.size(); }

    void
    addEulerianCell(InterfaceCellType* cell) {
      _supports.emplace_back(cell, 0.0);
    }

    Support&
    operator[](std::size_t const& index) { return _supports[index]; }

    Uni_PublicProperty(LagrangianNode*, lagrangianNode);
    Uni_PublicProperty(Supports,        supports);
  };

private:
  ScalarType
  computeImqRbf(VectorDsType one, VectorDsType two) const {
    auto radius = (one - two).norm();

    return computeImqRbf(radius);
  }

  ScalarType
  computeImqRbf(ScalarType radius) const {
    return 1.0 / std::sqrt(1.0 + _imqShape * _imqShape * radius * radius);
  }

  precice::SolverInterface* _preciceInterface;
  MemoryType const*         _memory;

  std::set<int> _bodyMeshSet;

  bool       _doComputeSupportRadius;
  ScalarType _suppportRadius;
  bool       _doComputeImqShape;
  ScalarType _imqShape;

  bool _isDevelopingStructure;

  Eigen::Matrix<unsigned, Dimensions, 1> _localEulerianCellSize;
  Eigen::Matrix<unsigned, Dimensions, 1> _globalEulerianCellSize;

  Set _velocitySet;
  Set _forceSet;

  int                         _displacementsId;
  std::vector<LagrangianNode> _lagrangianNodes;
  std::vector<int>            _lagrangianIds;
  std::vector<ScalarType>     _lagrangianDisplacements;
  std::array<std::vector<LagrangianNodeSupports>,
             Dimensions> _lagrangianNodesSupports;

  unsigned                                      _lagrangianIdOffset;
  std::array<unsigned, Dimensions>              _eulerianIdOffset;
  std::array<std::map<int, Set>, Dimensions>    _foreignEulerianCells;
  std::array<std::vector<unsigned>, Dimensions> receive_sizes;
  std::array < std::vector < std::vector<unsigned >>, Dimensions>
  receive_serial_indices;
  std::array<unsigned, Dimensions> send_size;
  std::array<unsigned, Dimensions> receive_size;
};
}
}
}
