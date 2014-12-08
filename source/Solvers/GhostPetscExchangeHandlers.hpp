#ifndef FsiSimulation_Ghost_PetscExchange_Handler_hpp
#define FsiSimulation_Ghost_PetscExchange_Handler_hpp

#include "GhostHandlersUtilities.hpp"
#include "StructuredMemory/Pointers.hpp"

#include <Uni/Logging/macros>

#include <petscdm.h>
#include <petscdmda.h>

#include <functional>
#include <memory>

namespace FsiSimulation {
namespace Solvers {
namespace Ghost {
namespace PetscExchange {
template <int D>
using Functor =
        std::function
        <void(typename StructuredMemory::Pointers<PetscScalar, D>::Type array)>;
template <int D>
using FunctorStack = FunctorStack<Functor<D>, D>;

template <int D>
inline Functor<D>
getEmptyFunctor() {
  return Functor<D>([] (
                      typename StructuredMemory::Pointers<PetscScalar, D>::Type
                      array) {});
}

template <typename TGrid,
          typename TAction,
          int D,
          int TDimension,
          int TDirection>
class Handler {
private:
  typedef StructuredMemory::Pointers<PetscScalar, D> Pointers;

public:
  Handler(TGrid const* grid,
          TAction*     action) : _grid(grid),
                                 _action(action) {}

  ~Handler() {}

  static Functor<D>
  getHandler(TGrid const* grid,
             TAction*     action) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, action));

    return Functor<D>(std::bind(&Handler::initialize, pointer, _1));
  }

  void
  initialize(typename Pointers::Type array) {
    for (auto const& accessor :
         _grid->indentedBoundaries[TDimension][TDirection]) {
      _action->exchange(accessor, array);
    }
  }

private:
  TGrid const* _grid;
  TAction*     _action;
};
}
}
}
}

#endif
