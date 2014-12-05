#ifndef FsiSimulation_Ghost_Velocity_Handler_hpp
#define FsiSimulation_Ghost_Velocity_Handler_hpp

#include "GhostHandlersUtilities.hpp"
#include "ParallelTopology.hpp"

#include <Uni/Logging/macros>

#include <functional>
#include <memory>

namespace FsiSimulation {
namespace Solvers {
namespace Ghost {
namespace Initialization {
typedef std::function<void ()> Functor;
template <int D>
using FunctorStack = FunctorStack<Functor, D>;

template <typename TGrid,
          typename Scalar,
          typename TAction,
          int D,
          int TDimension,
          int TDirection>
class Handler {
public:
  typedef ParallelTopology<D> SpecializedParallelTopology;

  Handler(TGrid const*                       grid,
          SpecializedParallelTopology const* parallelTopology,
          TAction*                           action)
    : _grid(grid),
      _parallelTopology(parallelTopology),
      _action(action) {}

  ~Handler() {}

  static Functor
  getHandler(TGrid const*                       grid,
             SpecializedParallelTopology const* parallelTopology,
             TAction*                           action) {
    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology, action));

    return Functor(std::bind(&Handler::initialize, pointer));
  }

  void
  initialize() {
    for (auto const& accessor : _grid->boundaries[TDimension][TDirection]) {
      auto newAccessor = accessor;
      newAccessor.initialize(accessor.relativeIndex(TDimension, !TDirection));
      //INFO << accessor.relativeIndex(TDimension, !TDirection).transpose();
      _action->setValue(accessor, newAccessor);
    }
  }

private:
  TGrid const*                       _grid;
  SpecializedParallelTopology const* _parallelTopology;
  TAction*                           _action;
};

template <typename TGrid,
          typename Scalar,
          typename TAction,
          int D>
class Stack {};

template <typename TGrid,
          typename Scalar,
          typename TAction>
class Stack<TGrid, Scalar, TAction, 2> {
public:
  typedef ParallelTopology<2> SpecializedParallelTopology;
  template <int TDimension, int TDirection>
  using _Handler = Handler<TGrid,
                           Scalar,
                           TAction,
                           2,
                           TDimension,
                           TDirection>;

  static FunctorStack<2>
  create(TGrid const*                       grid,
         SpecializedParallelTopology const* topology,
         TAction const*                     action) {
    FunctorStack<2> _instance = {
      _Handler<0, 0>::getHandler(grid, topology, action),
      _Handler<0, 1>::getHandler(grid, topology, action),
      _Handler<1, 0>::getHandler(grid, topology, action),
      _Handler<1, 1>::getHandler(grid, topology, action)
    };

    return _instance;
  }
};

template <typename TGrid,
          typename Scalar,
          typename TAction>
class Stack<TGrid, Scalar, TAction, 3> {
public:
  typedef ParallelTopology<3> SpecializedParallelTopology;
  template <int TDimension, int TDirection>
  using _Handler = Handler<TGrid,
                           Scalar,
                           TAction,
                           3,
                           TDimension,
                           TDirection>;

  static FunctorStack<3>
  create(TGrid const*                       grid,
         SpecializedParallelTopology const* topology,
         TAction const*                     action) {
    FunctorStack<3> _instance = {
      _Handler<0, 0>::getHandler(grid, topology, action),
      _Handler<0, 1>::getHandler(grid, topology, action),
      _Handler<1, 0>::getHandler(grid, topology, action),
      _Handler<1, 1>::getHandler(grid, topology, action),
      _Handler<2, 0>::getHandler(grid, topology, action),
      _Handler<2, 1>::getHandler(grid, topology, action)
    };

    return _instance;
  }
};
}
}
}
}

#endif
