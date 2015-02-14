#ifndef FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Actions_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_PetscExchange_Actions_hpp

#include "Private/utilities.hpp"

#include "FluidSimulation/Configuration.hpp"
#include "FluidSimulation/ParallelDistribution.hpp"

#include "StructuredMemory/Pointers.hpp"

#include <Uni/Logging/macros>

#include <petscdm.h>
#include <petscdmda.h>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace PetscExchange {
template <int TSolverDimension,
          int TDimension,
          int TDirection>
class VpeInputRhsGenerationAction {
public:
  VpeInputRhsGenerationAction(Configuration* configuration)
    : _configuration(configuration) {}

  template <typename TCellAccessor>
  inline void
  exchange(
    typename StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::CellType::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    typedef TCellAccessor                           CellAccessorType;
    typedef typename CellAccessorType::CellType     CellType;
    typedef typename CellAccessorType::VectorDiType VectorDiType;
    typedef typename CellType::Scalar               Scalar;

    enum {
      Dimensions = CellType::Dimensions
    };

    typedef StructuredMemory::Pointers<PetscScalar, Dimensions> Pointers;
    Pointers::dereference(array, index)
      = _configuration->walls[TDimension][TDirection]
        ->velocity() (TSolverDimension);
  }

private:
  Configuration* _configuration;
};

template <int TSolverDimension,
          int TDimension,
          int TDirection>
class VpeParabolicInputRhsGenerationAction {
public:
  VpeParabolicInputRhsGenerationAction(Configuration* configuration)
    : _configuration(configuration) {}

  template <typename TCellAccessor>
  inline void
  exchange(
    typename StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::CellType::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    typedef TCellAccessor                           CellAccessorType;
    typedef typename CellAccessorType::CellType     CellType;
    typedef typename CellAccessorType::VectorDiType VectorDiType;
    typedef typename CellType::Scalar               Scalar;

    enum {
      Dimensions = CellType::Dimensions
    };

    typedef StructuredMemory::Pointers<PetscScalar, Dimensions> Pointers;

    Scalar temp;

    if (TSolverDimension == TDimension) {
      temp = computeParabolicInputVelocity(
        accessor,
        _configuration->walls[TDimension][TDirection]
        ->velocity(),
        TSolverDimension);
    } else {
      temp = _configuration->walls[TDimension][TDirection]
             ->velocity() (TSolverDimension);
    }

    Pointers::dereference(array, index) = temp;

  }

private:
  Configuration* _configuration;
};

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
           VectorDiType const&     index,
           CellAccessorType const& accessor) {
    Pointers::dereference(array, index) = _value;
  }

private:
  Scalar _value;
};

template <int TSolverDimension>
class VpeRhsAcquiererAction {
public:
  template <typename TCellAccessor>
  inline void
  exchange(
    typename StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::CellType::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    typedef TCellAccessor                       CellAccessorType;
    typedef typename CellAccessorType::CellType CellType;
    typedef StructuredMemory::Pointers
      <PetscScalar, CellType::Dimensions> Pointers;
    accessor.currentCell()->velocity(TSolverDimension)
      = Pointers::dereference(array, index);
  }
};

class PpeRhsAcquiererAction {
public:
  template <typename TCellAccessor>
  inline void
  exchange(
    typename StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::CellType::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
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
