#ifndef FsiSimulation_FluidSimulation_GhostLayer_MpiExchange_Handler_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_MpiExchange_Handler_hpp

#include "Private/utilities.hpp"

#include "FluidSimulation/ParallelDistribution.hpp"
#include "FluidSimulation/Private/mpigenerics.hpp"

#include "StructuredMemory/Memory.hpp"

#include <Uni/Logging/macros>

#include <functional>
#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace MpiExchange {
typedef std::function<void (MPI_Comm&)> Functor;
template <int TD>
using FunctorStack = FunctorStack<Functor, TD>;

template <int TD>
Functor
getEmptyFunctor() {
  return Functor([] (MPI_Comm&) {});
}

template <typename TCell, typename TGrid>
using GetValueFunction = TCell * (*)(typename TGrid::CellAccessor const &);

template <typename TCell,
          typename TGrid,
          GetValueFunction<TCell, TGrid> TgetValue,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class Handler {
public:
  typedef ParallelDistribution<TD>                SpecializedParallelTopology;
  typedef StructuredMemory::Memory<TCell, TD - 1> Memory;
  typedef typename Memory::VectorDi               InternalVectorDi;
  typedef typename TGrid::VectorDi                VectorDi;

  Handler(TGrid const*                       grid,
          SpecializedParallelTopology const* parallelTopology,
          VectorDi const&                    leftIndent = VectorDi::Zero(),
          VectorDi const&                    rightIndent = VectorDi::Zero())
    : _ghostGrid(*grid),
      _boundaryGrid(*grid),
      _parallelTopology(parallelTopology) {
    InternalVectorDi size;
    VectorDi         ghostGridLeftIndent(VectorDi::Zero());
    VectorDi         ghostGridRightIndent(VectorDi::Zero());
    VectorDi         boundaryGridLeftIndent(VectorDi::Zero());
    VectorDi         boundaryGridRightIndent(VectorDi::Zero());
    int              i = 0;

    for (int d = 0; d < TD; ++d) {
      if (d != TDimension) {
        size(i) = grid->size(d) - leftIndent(d) - rightIndent(d);

        ghostGridLeftIndent(d)     = leftIndent(d);
        ghostGridRightIndent(d)    = rightIndent(d);
        boundaryGridLeftIndent(d)  = rightIndent(d);
        boundaryGridRightIndent(d) = leftIndent(d);

        ++i;
      }
    }
    _memory.allocate(size);

    if (TDirection == 0) {
      ghostGridLeftIndent(TDimension)  = 0;
      ghostGridRightIndent(TDimension) =
        grid->size(TDimension) - leftIndent(TDimension);
      boundaryGridLeftIndent(TDimension)  = leftIndent(TDimension);
      boundaryGridRightIndent(TDimension) =
        grid->size(TDimension) - 2 * leftIndent(TDimension);
    } else {
      ghostGridLeftIndent(TDimension) =
        grid->size(TDimension) - rightIndent(TDimension);
      ghostGridRightIndent(TDimension)   = 0;
      boundaryGridLeftIndent(TDimension) =
        grid->size(TDimension) - 2 * rightIndent(TDimension);
      boundaryGridRightIndent(TDimension) = rightIndent(TDimension);
    }
    _ghostGrid.setIndents(ghostGridLeftIndent, ghostGridRightIndent);
    _boundaryGrid.setIndents(boundaryGridLeftIndent, boundaryGridRightIndent);
  }

  ~Handler() {}

  static Functor
  getSendHandler(
    TGrid const*                       grid,
    SpecializedParallelTopology const* parallelTopology,
    VectorDi const&                    leftIndent = VectorDi::Zero(),
    VectorDi const&                    rightIndent = VectorDi::Zero()) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid,
                                  parallelTopology,
                                  leftIndent,
                                  rightIndent));

    return Functor(std::bind(&Handler::send, pointer, _1));
  }

  static Functor
  getReceiveHandler(
    TGrid const*                       grid,
    SpecializedParallelTopology const* parallelTopology,
    VectorDi const&                    leftIndent = VectorDi::Zero(),
    VectorDi const&                    rightIndent = VectorDi::Zero()) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid,
                                  parallelTopology,
                                  leftIndent,
                                  rightIndent));

    return Functor(std::bind(&Handler::receive, pointer, _1));
  }

  static Functor
  getExchangeHandler(
    TGrid const*                       grid,
    SpecializedParallelTopology const* parallelTopology,
    VectorDi const&                    leftIndent = VectorDi::Zero(),
    VectorDi const&                    rightIndent = VectorDi::Zero()) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid,
                                  parallelTopology,
                                  leftIndent,
                                  rightIndent));

    return Functor(std::bind(&Handler::exchange, pointer, _1));
  }

  void
  send(MPI_Comm& mpiCommunicator) {
    for (auto const& accessor : _boundaryGrid) {
      InternalVectorDi index;
      int              i = 0;

      for (int d = 0; d < TD; ++d) {
        if (d != TDimension) {
          index(i) = accessor.indexValue(d) - _boundaryGrid.leftIndent(d);
          ++i;
        }
      }
      *_memory.getCell(index) = *TgetValue(accessor);
    }
    MPI_Isend(
      _memory.data(),
      _memory.size().prod() * sizeof (TCell) / sizeof (TScalar),
      getMpiScalarType<TScalar>(),
      _parallelTopology->neighbors[TDimension][TDirection],
      0,
      mpiCommunicator,
      &_mpiRequest);
  }

  void
  receive(MPI_Comm& mpiCommunicator) {
    MPI_Recv(
      _memory.data(),
      _memory.size().prod() * sizeof (TCell) / sizeof (TScalar),
      getMpiScalarType<TScalar>(),
      _parallelTopology->neighbors[TDimension][TDirection],
      0,
      mpiCommunicator,
      &_mpiStatus);

    for (auto const& accessor : _ghostGrid) {
      InternalVectorDi index;
      int              i = 0;

      for (int d = 0; d < TD; ++d) {
        if (d != TDimension) {
          index(i) = accessor.indexValue(d) - _ghostGrid.leftIndent(d);
          ++i;
        }
      }
      *TgetValue(accessor) = *_memory.getCell(index);
    }
  }

  void
  exchange(MPI_Comm& mpiCommunicator) {
    for (auto const& accessor : _boundaryGrid) {
      InternalVectorDi index;
      int              i = 0;

      for (int d = 0; d < TD; ++d) {
        if (d != TDimension) {
          index(i) = accessor.indexValue(d) - _boundaryGrid.leftIndent(d);
          ++i;
        }
      }
      *_memory.getCell(index) = *TgetValue(accessor);
    }
    MPI_Sendrecv_replace(
      _memory.data(),
      _memory.size().prod() * sizeof (TCell) / sizeof (TScalar),
      getMpiScalarType<TScalar>(),
      _parallelTopology->neighbors[TDimension][TDirection],
      0,
      _parallelTopology->neighbors[TDimension][TDirection],
      0,
      mpiCommunicator,
      &_mpiStatus);

    for (auto const& accessor : _ghostGrid) {
      InternalVectorDi index;
      int              i = 0;

      for (int d = 0; d < TD; ++d) {
        if (d != TDimension) {
          index(i) = accessor.indexValue(d) - _ghostGrid.leftIndent(d);
          ++i;
        }
      }
      *TgetValue(accessor) = *_memory.getCell(index);
    }
  }

public:
  template <typename TType>
  using TwoDArray = std::array<std::array<TType, 2>, TD>;

private:
  TGrid                              _ghostGrid;
  TGrid                              _boundaryGrid;
  SpecializedParallelTopology const* _parallelTopology;
  Memory                             _memory;
  MPI_Request                        _mpiRequest;
  MPI_Status                         _mpiStatus;
};
}
}
}
}
#endif
