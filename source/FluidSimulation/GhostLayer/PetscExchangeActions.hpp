#ifndef FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Actions_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Actions_hpp

#include "FluidSimulation/ParallelDistribution.hpp"

#include "StructuredMemory/Pointers.hpp"

#include <petscdm.h>
#include <petscdmda.h>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace PetscExchange {
template <typename TCellAccessor>
class ConstantRhsGenerationAction {
public:
  typedef TCellAccessor                           CellAccessorType;
  typedef typename CellAccessorType::CellType     CellType;
  typedef typename CellAccessorType::VectorDiType VectorDiType;
  typedef typename CellType::Scalar               Scalar;

  enum {
    Dimensions = CellType::Dimensions
  };

private:
  typedef StructuredMemory::Pointers<PetscScalar, Dimensions> Pointers;

public:
  ConstantRhsGenerationAction(Scalar const& value)
    : _value(value) {}

  inline void
  exchange(typename Pointers::Type array,
           VectorDiType const&              index,
           CellAccessorType const& accessor) {
    Pointers::dereference(array, index) = _value;
  }

private:
  Scalar _value;
};

class PpeRhsAcquiererAction {
public:
  template <typename TCellAccessor>
  inline void
  exchange(
    typename StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::CellType::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const&
    accessor) {
    typedef TCellAccessor                       CellAccessorType;
    typedef typename CellAccessorType::CellType CellType;
    typedef StructuredMemory::Pointers
      <PetscScalar, CellType::Dimensions> Pointers;
    accessor.currentCell()->pressure() = Pointers::dereference(array, index);
  }
};
}
}
}
}

#endif
