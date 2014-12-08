#ifndef FsiSimulation_Ghost_PetscExchange_Actions_hpp
#define FsiSimulation_Ghost_PetscExchange_Actions_hpp

#include "ParallelTopology.hpp"
#include "Parameters.h"
#include "StructuredMemory/Pointers.hpp"
#include <petscdm.h>
#include <petscdmda.h>

namespace FsiSimulation {
namespace Solvers {
namespace Ghost {
namespace PetscExchange {
template <typename TGrid,
          int D,
          int TDimension,
          int TDirection>
class MovingWallRhsAction {
public:
  typedef StructuredMemory::Pointers<PetscScalar, D> Pointers;
  typedef ParallelTopology<D>                        SpecializedParallelTopology;

public:
  MovingWallRhsAction(
    Parameters&                        parameters,
    SpecializedParallelTopology const* parallelTopology)
    : _parameters(parameters),
      _parallelTopology(parallelTopology) {}

  void
  exchange(typename TGrid::CellAccessor const& current,
           typename Pointers::Type             array) {
    auto corner = _parallelTopology->corner;

    auto index = current.indexValues();
    index += corner;

    Pointers::dereference(array, index) = 0.0;
  }

  Parameters&                        _parameters;
  SpecializedParallelTopology const* _parallelTopology;
};

template <typename TGrid,
          int D,
          int TDimension,
          int TDirection>
class CopyPressureAction {
public:
  typedef StructuredMemory::Pointers<PetscScalar, D> Pointers;
  typedef ParallelTopology<D>                        SpecializedParallelTopology;

public:
  CopyPressureAction(SpecializedParallelTopology const* parallelTopology)
    : _parallelTopology(parallelTopology) {}

  void
  exchange(typename TGrid::CellAccessor const& current,
           typename Pointers::Type             array) {
    auto corner = _parallelTopology->corner;

    auto index = current.indexValues();
    index += corner;

    current.currentCell()->pressure() = Pointers::dereference(array, index);
  }

  SpecializedParallelTopology const* _parallelTopology;
};
}
}
}
}

#endif
