#pragma once

#include "Private/utilities.hpp"

#include "Simulation/ParallelDistribution.hpp"
#include "Simulation/Private/mpigenerics.hpp"

#include <Uni/Logging/macros>

#include <functional>
#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace MpiExchange {
typedef std::function<void ()> Functor;
template <int Dimensions>
using FunctorStack = BasicFunctorStack<Functor, Dimensions>;

template <int Dimensions>
Functor
getEmptyFunctor() {
  return Functor([] () {});
}

template <typename TGrid,
          typename TDataElement>
using GetValueFunction
        = TDataElement * (*)(typename TGrid::CellAccessor const &);

template <typename TGrid,
          typename TDataElement>
using SetValueFunction
        = void (*)(typename TGrid::CellAccessor const&,
                   int const&,
                   TDataElement const&);

template <typename TDataElement,
          unsigned TDataElementSize,
          typename TGrid,
          GetValueFunction<TGrid, TDataElement> TgetValue,
          SetValueFunction<TGrid, TDataElement> TsetValue,
          int TDimension,
          int TDirection>
class Handler {
public:
  using GridType = TGrid;

  using CellAccessorType = typename GridType::CellAccessor;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

  using ScalarType = typename CellAccessorType::ScalarType;

  using ParallelDistributionType =  ParallelDistribution<Dimensions>;

  using InternalVectorDi = Eigen::Matrix<int, Dimensions - 1, 1>;

  using VectorDiType = typename TGrid::VectorDi;

  Handler(TGrid const*                    grid,
          ParallelDistributionType const* parallelTopology)
    : _ghostGrid(*grid),
    _boundaryGrid(*grid),
    _parallelDistribution(parallelTopology) {
    InternalVectorDi size;
    VectorDiType     ghostGridLeftIndent(VectorDiType::Zero());
    VectorDiType     ghostGridRightIndent(VectorDiType::Zero());
    VectorDiType     boundaryGridLeftIndent(VectorDiType::Zero());
    VectorDiType     boundaryGridRightIndent(VectorDiType::Zero());
    int              i = 0;

    for (int d = 0; d < Dimensions; ++d) {
      if (d != TDimension) {
        size(i) = _ghostGrid.innerSize(d);

        ghostGridLeftIndent(d)     = _ghostGrid.leftIndent(d);
        ghostGridRightIndent(d)    = _ghostGrid.rightIndent(d);
        boundaryGridLeftIndent(d)  = _ghostGrid.leftIndent(d);
        boundaryGridRightIndent(d) = _ghostGrid.rightIndent(d);

        ++i;
      }
    }

    _rowMemorySize = size.prod() * TDataElementSize
                     * _ghostGrid.indent(TDirection)(TDimension);
    _rowMemory.reset(new ScalarType[_rowMemorySize]);

    if (TDirection == 0) {
      ghostGridLeftIndent(TDimension) = 0;
      ghostGridRightIndent(TDimension)
        = _ghostGrid.size(TDimension) -
          _ghostGrid.leftIndent(TDimension);
      boundaryGridLeftIndent(TDimension)
        = _boundaryGrid.leftIndent(TDimension);
      boundaryGridRightIndent(TDimension)
        = _boundaryGrid.size(TDimension)
          - 2 * _ghostGrid.leftIndent(TDimension);
    } else {
      ghostGridLeftIndent(TDimension)
        = _ghostGrid.size(TDimension)
          - _ghostGrid.rightIndent(TDimension);
      ghostGridRightIndent(TDimension) = 0;
      boundaryGridLeftIndent(TDimension)
        = _boundaryGrid.size(TDimension)
          - 2 * _ghostGrid.rightIndent(TDimension);
      boundaryGridRightIndent(TDimension)
        = _boundaryGrid.rightIndent(TDimension);
    }
    _ghostGrid.setIndents(ghostGridLeftIndent, ghostGridRightIndent);
    _boundaryGrid.setIndents(boundaryGridLeftIndent, boundaryGridRightIndent);
  }

  ~Handler() {}

  static Functor
  getSendHandler(
    TGrid const*                    grid,
    ParallelDistributionType const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid,
                                  parallelTopology));

    return Functor(std::bind(&Handler::send, pointer));
  }

  static Functor
  getReceiveHandler(
    TGrid const*                    grid,
    ParallelDistributionType const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid,
                                  parallelTopology));

    return Functor(std::bind(&Handler::receive, pointer));
  }

  static Functor
  getExchangeHandler(
    TGrid const*                    grid,
    ParallelDistributionType const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid,
                                  parallelTopology));

    return Functor(std::bind(&Handler::exchange, pointer));
  }

  void
  send() {
    int index = 0;

    for (auto const& accessor : _boundaryGrid) {
      for (unsigned i = 0; i < TDataElementSize; ++i) {
        _rowMemory.get()[index] = TgetValue(accessor)[i];
        ++index;
      }
    }

    MPI_Isend(
      _rowMemory.get(),
      _rowMemorySize,
      Private::getMpiScalarType<TDataElement>(),
      _parallelDistribution->neighbors[TDimension][TDirection],
      1001,
      _parallelDistribution->mpiCommunicator,
      &_mpiRequest);
  }

  void
  receive() {
    MPI_Recv(
      _rowMemory.get(),
      _rowMemorySize,
      Private::getMpiScalarType<TDataElement>(),
      _parallelDistribution->neighbors[TDimension][TDirection],
      1001,
      _parallelDistribution->mpiCommunicator,
      &_mpiStatus);

    int index = 0;

    for (auto const& accessor : _ghostGrid) {
      for (unsigned i = 0; i < TDataElementSize; ++i) {
        TsetValue(accessor, i, _rowMemory.get()[index]);
        ++index;
      }
    }
  }

  void
  exchange() {
    int index = 0;

    for (auto const& accessor : _boundaryGrid) {
      for (unsigned i = 0; i < TDataElementSize; ++i) {
        _rowMemory.get()[index] = TgetValue(accessor)[i];
        ++index;
      }
    }

    MPI_Sendrecv_replace(
      _rowMemory.get(),
      _rowMemorySize,
      Private::getMpiScalarType<TDataElement>(),
      _parallelDistribution->neighbors[TDimension][TDirection],
      1002,
      _parallelDistribution->neighbors[TDimension][TDirection],
      1002,
      _parallelDistribution->mpiCommunicator,
      &_mpiStatus);

    index = 0;

    for (auto const& accessor : _ghostGrid) {
      for (unsigned i = 0; i < TDataElementSize; ++i) {
        TsetValue(accessor, i, _rowMemory.get()[index]);
        ++index;
      }
    }
  }

private:
  TGrid                           _ghostGrid;
  TGrid                           _boundaryGrid;
  ParallelDistributionType const* _parallelDistribution;
  std::unique_ptr<ScalarType[]>   _rowMemory;
  unsigned long long              _rowMemorySize;
  MPI_Request                     _mpiRequest;
  MPI_Status                      _mpiStatus;
};
}
}
}
}
