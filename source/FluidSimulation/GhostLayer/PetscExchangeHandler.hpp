#ifndef FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Handler_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Handler_hpp

#include "FluidSimulation/ParallelDistribution.hpp"

#include "Private/utilities.hpp"

#include "StructuredMemory/Pointers.hpp"

#include <Uni/Logging/macros>

#include <petscdm.h>
#include <petscdmda.h>

#include <functional>
#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace PetscExchange {
template <int TD>
using Functor
        = std::function
          <void(typename StructuredMemory::Pointers<PetscScalar, TD>::Type)>;
template <int TD>
using FunctorStack = FunctorStack<Functor<TD>, TD>;

template <int TD>
inline Functor<TD>
getEmptyFunctor() {
  return Functor<TD>(
    [] (typename StructuredMemory::Pointers<PetscScalar, TD>::Type array) {});
}

template <typename TGrid,
          typename TAction,
          int TDimension,
          int TDirection>
class Handler {
public:
  typedef TGrid                               GridType;
  typedef typename GridType::CellAccessorType CellAccessorType;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

  typedef ParallelDistribution<Dimensions> ParallelDistributionType;

private:
  typedef StructuredMemory::Pointers<PetscScalar, Dimensions> Pointers;

public:
  Handler(GridType const*                 grid,
          ParallelDistributionType const* parallelDistribution,
          TAction*                        action)
    : _grid(grid),
    _parallelDistribution(parallelDistribution),
    _action(action) {}

  ~Handler() {}

  static Functor<Dimensions>
  getHandler(GridType const*                 grid,
             ParallelDistributionType const* parallelDistribution,
             TAction*                        action) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelDistribution, action));

    return Functor<Dimensions>(std::bind(&Handler::initialize, pointer, _1));
  }

  void
  initialize(typename Pointers::Type array) {
    for (auto& accessor :
         _grid->indentedBoundaries[TDimension][TDirection]) {
      auto corner = _parallelDistribution->corner;

      auto index = accessor.indexValues();
      index += corner;

      _action->exchange(array, index, accessor);
    }
  }

private:
  TGrid const*                    _grid;
  ParallelDistributionType const* _parallelDistribution;
  TAction*                        _action;
};
}
}
}
}

#endif
