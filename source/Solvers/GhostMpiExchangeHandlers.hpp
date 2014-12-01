#ifndef FsiSimulation_Solvers_GhostMpiExchangeHandlers_hpp
#define FsiSimulation_Solvers_GhostMpiExchangeHandlers_hpp

#include "GhostHandlersUtilities.hpp"
#include "ParallelTopology.hpp"
#include "StructuredMemory/Memory.hpp"
#include "mpigenerics.hpp"

#include <Uni/Logging/macros>

#include <functional>
#include <memory>

namespace FsiSimulation {
namespace Solvers {
namespace Ghost {
namespace MpiExchange {
typedef std::function<void (MPI_Comm&)> GetFalueFunctor;
typedef std::function<void (MPI_Comm&)> Functor;
template <int D>
using FunctorStack = FunctorStack<Functor, D>;

template <typename TCell, typename TGrid>
using GetValueFunction = TCell * (*)(typename TGrid::CellAccessor const &);

template <typename TCell,
          typename TGrid,
          typename Scalar,
          int D,
          int TDimension,
          int TDirection,
          GetValueFunction<TCell, TGrid> TgetValue>
class Handler {
public:
  typedef ParallelTopology<D>                    SpecializedParallelTopology;
  typedef StructuredMemory::Memory<TCell, D - 1> Memory;
  typedef typename Memory::VectorDi              VectorDi;

  Handler(TGrid const*                       grid,
          SpecializedParallelTopology const* parallelTopology)
    : _grid(grid),
      _parallelTopology(parallelTopology) {
    VectorDi indexShift;
    VectorDi size;
    int      i = 0;

    for (int d = 0; d < D; ++d) {
      if (d != TDimension) {
        indexShift(i) = grid->indent (TDirection)(d);
        size(i)       = grid->size(d) - indexShift(i);
      }
    }
    _memory.indexShift(indexShift);
    _memory.allocate(size);
  }

  ~Handler() {
    logInfo("Mpi Exchange Handler is Destroyed");
  }

  static Functor
  getHandler(TGrid const*                       grid,
             SpecializedParallelTopology const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology));

    return Functor(std::bind(&Handler::exchange,
                             pointer,
                             _1));
  }

  void
  exchange(MPI_Comm& mpiCommunicator) {
    for (auto const& accessor : _grid->boundaries[TDimension][TDirection]) {
      VectorDi index;
      int      i = 0;

      for (int d = 0; d < D; ++d) {
        if (d != TDimension) {
          index(i) = accessor.indexValue(d);
          ++i;
        }
      }
      *_memory.getCell(index) = *TgetValue(accessor);
    }
    MPI_Sendrecv_replace(
      _memory.data(),
      _memory.size().prod(),
      getMpiScalarType<Scalar>(),
      _parallelTopology->neighbors[TDimension][TDirection],
      0,
      _parallelTopology->neighbors[TDimension][!TDirection],
      0,
      mpiCommunicator,
      &_mpiStatus);

    for (auto const& accessor : _grid->boundaries[TDimension][!TDirection]) {
      VectorDi index;
      int      i = 0;

      for (int d = 0; d < D; ++d) {
        if (d != TDimension) {
          index(i) = accessor.indexValue(d);
          ++i;
        }
      }
      *TgetValue(accessor) = *_memory.getCell(index);
    }
  }

public:
  template <typename TType>
  using TwoDArray = std::array<std::array<TType, 2>, D>;

private:
  TGrid const*                       _grid;
  SpecializedParallelTopology const* _parallelTopology;
  Memory                             _memory;
  MPI_Status                         _mpiStatus;
};

template <typename TCell,
          typename TGrid,
          typename Scalar,
          int D,
          TCell* (* getValue)(typename TGrid::CellAccessor const&)>
class Stack {};

template <typename TCell,
          typename TGrid,
          typename Scalar,
          TCell* (* getValue)(typename TGrid::CellAccessor const&)>
class Stack<TCell, TGrid, Scalar, 2, getValue> {
public:
  typedef ParallelTopology<2> SpecializedParallelTopology;
  template <int TDimension, int TDirection>
  using _Handler = Handler<TCell,
                           TGrid,
                           Scalar,
                           2,
                           TDimension,
                           TDirection,
                           getValue>;

  static FunctorStack<2>
  create(TGrid const*                       grid,
         SpecializedParallelTopology const* topology) {
    FunctorStack<2> _instance =
    { _Handler<0, 0>::getHandler(grid, topology),
      _Handler<0, 1>::getHandler(grid, topology),
      _Handler<1, 0>::getHandler(grid, topology),
      _Handler<1, 1>::getHandler(grid, topology) };

    return _instance;
  }
};

template <typename TCell,
          typename TGrid,
          typename Scalar,
          TCell* (* getValue)(typename TGrid::CellAccessor const&)>
class Stack<TCell, TGrid, Scalar, 3, getValue> {
public:
  typedef ParallelTopology<3> SpecializedParallelTopology;
  template <int TDimension, int TDirection>
  using _Handler = Handler<TCell,
                           TGrid,
                           Scalar,
                           3,
                           TDimension,
                           TDirection,
                           getValue>;

  static FunctorStack<3>
  create(TGrid const*                       grid,
         SpecializedParallelTopology const* topology) {
    FunctorStack<3> _instance =
    { _Handler<0, 0>::getHandler(grid, topology),
      _Handler<0, 1>::getHandler(grid, topology),
      _Handler<1, 0>::getHandler(grid, topology),
      _Handler<1, 1>::getHandler(grid, topology),
      _Handler<2, 0>::getHandler(grid, topology),
      _Handler<2, 1>::getHandler(grid, topology) };

    return _instance;
  }
};
}
}
}
}
#endif
