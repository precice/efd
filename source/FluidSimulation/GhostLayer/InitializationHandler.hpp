#ifndef FsiSimulation_FluidSimulation_GhostLayer_Initialization_Handler_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_Initialization_Handler_hpp

#include "Private/utilities.hpp"

#include "FluidSimulation/ParallelDistribution.hpp"

#include <Uni/Logging/macros>

#include <functional>
#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace Initialization {
typedef std::function<void ()> Functor;
template <int TD>
using FunctorStack = FunctorStack<Functor, TD>;

template <int TD>
inline Functor
getEmptyFunctor() {
  return Functor([] () {});
}

template <typename TGrid,
          typename TAction,
          int TD,
          int TDimension,
          int TDirection>
class Handler {
public:
  typedef ParallelDistribution<TD> SpecializedParallelTopology;

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
    for (auto const& accessor :* _grid) {
      auto newAccessor = accessor;
      newAccessor.initialize(accessor.relativeIndex(TDimension, !TDirection));
      _action->setValue(accessor, newAccessor);
    }
  }

private:
  TGrid const*                       _grid;
  SpecializedParallelTopology const* _parallelTopology;
  TAction*                           _action;
};
}
}
}
}

#endif
