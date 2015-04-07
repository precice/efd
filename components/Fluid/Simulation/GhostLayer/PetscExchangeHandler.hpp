#pragma once

#include "Simulation/ParallelDistribution.hpp"

#include "Private/utilities.hpp"

#include <Uni/Logging/macros>
#include <Uni/StructuredMemory/Pointers>

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
          <void(typename Uni::StructuredMemory::Pointers<PetscScalar, TD>::Type)>;
template <int TD>
using FunctorStack = FunctorStack<Functor<TD>, TD>;

template <int TD>
inline Functor<TD>
getEmptyFunctor() {
  return Functor<TD>(
    [] (typename Uni::StructuredMemory::Pointers<PetscScalar, TD>::Type) {});
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
  using Pointers = Uni::StructuredMemory::Pointers<PetscScalar, Dimensions>;

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
    for (auto& accessor : _grid->indentedBoundaries[TDimension][TDirection]) {
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
