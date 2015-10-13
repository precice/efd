#pragma once

#include "Private/utilities.hpp"

#include "Simulation/Private/mpigenerics.hpp"

#include <Uni/ExecutionControl/exception>
#include <Uni/Logging/macros>

#include <functional>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <typename TDataElement,
          unsigned TDataElementSize,
          typename TSolverTraits,
          int TDimension1,
          int TDirection1,
          int TDimension2,
          int TDirection2>
class CornerVelocityHandler {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using MemoryType = typename SolverTraitsType::MemoryType;

  using GridType = typename SolverTraitsType::GridType;

  using SubgridType = typename SolverTraitsType::SubgridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

public:
  CornerVelocityHandler(
    SubgridType const*              grid,
    ParallelDistributionType const* parallel_distribution)
    : _grid(grid),
    _parallelDistribution(parallel_distribution),
    _dimension(-1) {
    int size = 1;

    for (unsigned d = 0; d < Dimensions; ++d) {
      if (d != TDimension1 && d != TDimension2) {
        size      *= _grid->innerSize(d);
        _dimension = d;
        break;
      }
    }

    _rowMemory.resize(size * TDataElementSize);

    _ghostIndex    = _grid->leftIndent();
    _boundaryIndex = _grid->leftIndent();
    VectorDiType communication_index = _parallelDistribution->index;
    composeIndexAndOffset(TDimension1, TDirection1, communication_index);
    composeIndexAndOffset(TDimension2, TDirection2, communication_index);
    _communicationRank
      = _parallelDistribution->getMpiRankFromSubdomainSpatialIndex(communication_index);

    // if (_communicationRank >= 0) {
    //   logInfo("@@@@@@@@@@ {1} {2} {3} | {4} {5}",
    //         TDimension1,
    //         TDimension2,
    //         _dimension,
    //         parallel_distribution->rank(),
    //         _communicationRank);
    // }
  }

  void
  computeCornerVelocity() {
    if (_communicationRank < 0) {
      return;
    }
    // logInfo("Enter");

    auto accessor = *_grid->at(_boundaryIndex);
    int  index    = 0;

    if (_dimension == -1) {
      for (unsigned i = 0; i < TDataElementSize; ++i) {
        _rowMemory[index] = accessor.velocity(i);
        ++index;
      }
    } else {
      for (; accessor.indexValue(_dimension)
           < _grid->innerLimit(_dimension);
           ++accessor.index(_dimension)) {
        for (unsigned i = 0; i < TDataElementSize; ++i) {
          _rowMemory[index] = accessor.velocity(i);
          ++index;
        }
      }
    }

    MPI_Status mpi_status;
    MPI_Sendrecv_replace(
      _rowMemory.data(),
      _rowMemory.size(),
      Private::getMpiScalarType<TDataElement>(),
      _communicationRank,
      1003,
      _communicationRank,
      1003,
      _parallelDistribution->mpiCommunicator,
      &mpi_status);

    accessor = *_grid->at(_ghostIndex);
    index    = 0;

    if (_dimension == -1) {
      for (unsigned i = 0; i < TDataElementSize; ++i) {
        accessor.velocity(i) = _rowMemory[index];
        ++index;
      }
    } else {
      for (; accessor.indexValue(_dimension)
           < _grid->innerLimit(_dimension);
           ++accessor.index(_dimension)) {
        for (unsigned i = 0; i < TDataElementSize; ++i) {
          // logInfo("{1} {2} {3}", _parallelDistribution->rank(),
          //         accessor.indexValue(_dimension),
          //         _rowMemory[index]);
          accessor.velocity(i) = _rowMemory[index];
          ++index;
        }
      }
    }
    // logInfo("Exit");
  }

private:
  void
  composeIndexAndOffset(int const&    dimension,
                        int const&    direction,
                        VectorDiType& communication_index) {
    if (direction == 0) {
      _ghostIndex(dimension)          = 0;
      _boundaryIndex(dimension)       = _grid->leftIndent(dimension);
      communication_index(dimension) -= 1;
    } else if (direction == 1) {
      _ghostIndex(dimension) = _grid->size(dimension)
                               - _grid->rightIndent(dimension);
      _boundaryIndex(dimension)
        = _grid->size(dimension)
          - _grid->rightIndent(dimension)
          - _grid->leftIndent(dimension);
      communication_index(dimension) += 1;
    } else {
      throwException("Unknown direction value");
    }
  }

  SubgridType const*              _grid;
  ParallelDistributionType const* _parallelDistribution;
  int                             _dimension;
  int                             _communicationRank;
  std::vector<ScalarType>         _rowMemory;
  VectorDiType                    _ghostIndex;
  VectorDiType                    _boundaryIndex;
};

template <typename TSolverTraits,
          int TDimensions = TSolverTraits::Dimensions>
struct CornerVelocityHandlers {};

template <typename TSolverTraits>
struct CornerVelocityHandlers<TSolverTraits, 2> {
  using SolverTraitsType = TSolverTraits;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using SubgridType = typename SolverTraitsType::SubgridType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using ScalarType = typename SolverTraitsType::ScalarType;
  CornerVelocityHandlers(SubgridType const*              grid,
                         ParallelDistributionType const* parallel_distribution)
    : _rb(grid, parallel_distribution),
    _lt(grid, parallel_distribution) {}

  void
  execute() {
    _rb.computeCornerVelocity();
    _lt.computeCornerVelocity();
  }

private:
  CornerVelocityHandler<ScalarType, 2, TSolverTraits, 0, 1, 1, 0> _rb;
  CornerVelocityHandler<ScalarType, 2, TSolverTraits, 0, 0, 1, 1> _lt;
};

template <typename TSolverTraits>
struct CornerVelocityHandlers<TSolverTraits, 3> {
  using SolverTraitsType = TSolverTraits;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using SubgridType = typename SolverTraitsType::SubgridType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  CornerVelocityHandlers(SubgridType const*              grid,
                         ParallelDistributionType const* parallel_distribution)
    // :
    // _rightBottom(grid, parallel_distribution),
    // _topLeft(grid, parallel_distribution)
    // _rightBack(grid, parallel_distribution),
    // _topBack(grid, parallel_distribution),
    // _frontLeft(grid, parallel_distribution),
    // _frontBottom(grid, parallel_distribution) 
    {
      (void)grid;
      (void)parallel_distribution;
    }

  void
  execute() {
    // _rightBottom.computeCornerVelocity();
    // _topLeft.computeCornerVelocity();
    
    // _rightBack.computeCornerVelocity();
    // _frontLeft.computeCornerVelocity();

    // _topBack.computeCornerVelocity();
    // _frontBottom.computeCornerVelocity();
  }

private:
  // CornerVelocityHandler<ScalarType, 3, TSolverTraits, 0, 1, 1, 0>
  // _rightBottom;
  // CornerVelocityHandler<ScalarType, 3, TSolverTraits, 1, 1, 0, 0>
  // _topLeft;
  // CornerVelocityHandler<ScalarType, 3, TSolverTraits, 0, 1, 2, 0>
  // _rightBack;
  // CornerVelocityHandler<ScalarType, 3, TSolverTraits, 1, 1, 2, 0>
  // _topBack;
  // CornerVelocityHandler<ScalarType, 3, TSolverTraits, 2, 1, 0, 0>
  // _frontLeft;
  // CornerVelocityHandler<ScalarType, 3, TSolverTraits, 2, 1, 1, 0>
  // _frontBottom;
};
}
}
}
