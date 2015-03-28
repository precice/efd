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
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    typedef TCellAccessor                           CellAccessorType;
    typedef typename CellAccessorType::VectorDiType VectorDiType;
    typedef typename CellAccessorType::ScalarType   ScalarType;

    enum {
      Dimensions = CellAccessorType::Dimensions
    };

    typedef StructuredMemory::Pointers<PetscScalar, Dimensions> Pointers;

    auto tempIndex = index;

    if (TSolverDimension == TDirection) {
      if (TDirection == 1) {
        tempIndex(TSolverDimension) -= 1;
      }
    }

    Pointers::dereference(array, tempIndex)
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
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    typedef TCellAccessor                           CellAccessorType;
    typedef typename CellAccessorType::VectorDiType VectorDiType;
    typedef typename CellAccessorType::ScalarType   ScalarType;

    enum {
      Dimensions = CellAccessorType::Dimensions
    };

    typedef StructuredMemory::Pointers<PetscScalar, Dimensions> Pointers;

    ScalarType temp;
    auto       tempIndex = index;

    temp = _configuration->walls[TDimension][TDirection]
           ->velocity()(TSolverDimension);

    if (TSolverDimension == TDimension) {
      compute_parabolic_input_velocity(
        accessor,
        TSolverDimension,
        temp);

      if (TDirection == 1) {
        tempIndex(TSolverDimension) -= 1;
      }
    }

    Pointers::dereference(array, tempIndex) = temp;
  }

private:
  Configuration* _configuration;
};

template <typename TCellAccessor,
          int TDimension>
class ConstantRhsGenerationAction {
public:
  typedef TCellAccessor                           CellAccessorType;
  typedef typename CellAccessorType::VectorDiType VectorDiType;
  typedef typename CellAccessorType::ScalarType   ScalarType;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

private:
  typedef StructuredMemory::Pointers<PetscScalar, Dimensions> Pointers;

public:
  ConstantRhsGenerationAction(ScalarType const& value,
                              int const&        offset = 0)
    : _value(value),
    _offset(offset) {}

  inline void
  exchange(typename Pointers::Type array,
           VectorDiType const&     index,
           CellAccessorType const& accessor) {
    // auto temp = index;
    // temp(TDimension)                  -= _offset;
    Pointers::dereference(array, index) = _value;
  }

private:
  ScalarType _value;
  int const  _offset;
};

template <int TSolverDimension>
class VpeRhsAcquiererAction {
public:
  template <typename TCellAccessor>
  inline void
  exchange(
    typename StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    typedef TCellAccessor CellAccessorType;

    typedef StructuredMemory::Pointers
      <PetscScalar, CellAccessorType::Dimensions> Pointers;

    accessor.fgh(TSolverDimension)
      = Pointers::dereference(array, index);
  }
};

class PpeRhsAcquiererAction1 {
public:
  template <typename TCellAccessor>
  inline void
  exchange(
    typename StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    typedef TCellAccessor CellAccessorType;
    typedef StructuredMemory::Pointers
      <PetscScalar, CellAccessorType::Dimensions> Pointers;
    accessor.pressure() = Pointers::dereference(array, index);
  }
};

class PpeRhsAcquiererAction2 {
public:
  template <typename TCellAccessor>
  inline void
  exchange(
    typename StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    typedef TCellAccessor CellAccessorType;
    typedef StructuredMemory::Pointers
      <PetscScalar, CellAccessorType::Dimensions> Pointers;
    auto const& value =  Pointers::dereference(array, index);

    // accessor.pressure()       =  accessor.projectionTerm() + value;
    accessor.projectionTerm() = value;
  }
};
}
}
}
}

#endif
