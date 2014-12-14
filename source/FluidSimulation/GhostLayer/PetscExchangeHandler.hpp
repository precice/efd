#ifndef FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Handler_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Handler_hpp

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
using Functor =
        std::function
        <void(typename StructuredMemory::Pointers<PetscScalar, TD>::Type
              array)>;
template <int TD>
using FunctorStack = FunctorStack<Functor<TD>, TD>;

template <int TD>
inline Functor<TD>
getEmptyFunctor() {
  return Functor<TD>([] (
                       typename StructuredMemory::Pointers<PetscScalar,
                                                           TD>::Type
                       array) {});
}

template <typename TGrid,
          typename TAction,
          int TD,
          int TDimension,
          int TDirection>
class Handler {
private:
  typedef StructuredMemory::Pointers<PetscScalar, TD> Pointers;

public:
  Handler(TGrid const* grid,
          TAction*     action) : _grid(grid),
                                 _action(action) {}

  ~Handler() {}

  static Functor<TD>
  getHandler(TGrid const* grid,
             TAction*     action) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, action));

    return Functor<TD>(std::bind(&Handler::initialize, pointer, _1));
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
