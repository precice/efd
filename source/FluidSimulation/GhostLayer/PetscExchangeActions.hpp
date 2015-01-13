#ifndef FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Actions_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Actions_hpp

#include "FluidSimulation/Configuration.hpp"
#include "FluidSimulation/ParallelDistribution.hpp"

#include "StructuredMemory/Pointers.hpp"

#include <petscdm.h>
#include <petscdmda.h>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace PetscExchange {
template <typename TGrid,
          int TD,
          int TDimension,
          int TDirection>
class MovingWallRhsAction {
public:
  typedef StructuredMemory::Pointers<PetscScalar, TD> Pointers;
  typedef ParallelDistribution<TD>
    SpecializedParallelTopology;

public:
  MovingWallRhsAction(
    Configuration*                     parameters,
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

  Configuration*                     _parameters;
  SpecializedParallelTopology const* _parallelTopology;
};

template <typename TGrid,
          int TD,
          int TDimension,
          int TDirection>
class CopyPressureAction {
public:
  typedef StructuredMemory::Pointers<PetscScalar, TD> Pointers;
  typedef ParallelDistribution<TD>
    SpecializedParallelTopology;

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
