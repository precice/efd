#pragma once

#include "Private/utilities.hpp"

#include "Simulation/Configuration.hpp"
#include "Simulation/ParallelDistribution.hpp"

#include <Uni/Logging/macros>
#include <Uni/StructuredMemory/Pointers>

#include <petscdm.h>
#include <petscdmda.h>

namespace Fluid {
namespace Simulation {
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
    typename Uni::StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    ((void)accessor);
    using CellAccessorType = TCellAccessor;

    enum {
      Dimensions = CellAccessorType::Dimensions
    };

    using Pointers = Uni::StructuredMemory::Pointers<PetscScalar, Dimensions>;

    auto tempIndex = index;

    if (TSolverDimension == TDimension) {
      if (TDirection == 1) {
        tempIndex(TSolverDimension)        -= 1;
        Pointers::dereference(array, index) = 0.0;
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
    typename Uni::StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    using CellAccessorType = TCellAccessor;

    using ScalarType = typename CellAccessorType::ScalarType;

    enum {
      Dimensions = CellAccessorType::Dimensions
    };

    using Pointers = Uni::StructuredMemory::Pointers
                     <PetscScalar, CellAccessorType::Dimensions>;

    ScalarType temp;
    auto       tempIndex = index;

    temp = _configuration->walls[TDimension][TDirection]
           ->velocity()(TSolverDimension);

    if (TSolverDimension == TDimension) {
      if (TDirection == 1) {
        tempIndex(TSolverDimension)        -= 1;
        Pointers::dereference(array, index) = 0.0;
        compute_parabolic_input_velocity(
          accessor.relative(TSolverDimension, -1),
          TSolverDimension,
          temp);
      } else {
        compute_parabolic_input_velocity(
          accessor,
          TSolverDimension,
          temp);
      }
    }

    Pointers::dereference(array, tempIndex) = temp;
  }

private:
  Configuration* _configuration;
};

template <typename TScalar,
          int TSolverDimension,
          int TDimension,
          int TDirection>
class VpeConstantRhsGenerationAction {
public:
  VpeConstantRhsGenerationAction(TScalar const& value)
    : _value(value) {}

  template <typename TCellAccessor>
  inline void
  exchange(
    typename Uni::StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const&                        accessor) {
    ((void)accessor);
    using Pointers = Uni::StructuredMemory::Pointers
                     <PetscScalar, TCellAccessor::Dimensions>;

    auto temp = index;

    if (TDimension == TSolverDimension) {
      if (TDirection == 1) {
        temp(TSolverDimension)             -= 1;
        Pointers::dereference(array, index) = 0.0;
      }
    }

    Pointers::dereference(array, temp) = _value;
  }

private:
  TScalar _value;
};

template <typename TScalar>
class PpeConstantRhsGenerationAction {
public:
  PpeConstantRhsGenerationAction(TScalar const& value)
    : _value(value) {}

  template <typename TCellAccessor>
  inline void
  exchange(
    typename Uni::StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const&                        accessor) {
    ((void)accessor);
    using Pointers = Uni::StructuredMemory::Pointers
                     <PetscScalar, TCellAccessor::Dimensions>;

    Pointers::dereference(array, index) = _value;
  }

private:
  TScalar _value;
};

template <int TSolverDimension,
          int TDimension,
          int TDirection>
class VpeRhsAcquiererAction {
public:
  template <typename TCellAccessor>
  inline void
  exchange(
    typename Uni::StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    using  CellAccessorType = TCellAccessor;

    using Pointers = Uni::StructuredMemory::Pointers
                     <PetscScalar, CellAccessorType::Dimensions>;

    accessor.fgh(TSolverDimension)
      = Pointers::dereference(array, index);
  }
};

class PpeRhsAcquiererAction1 {
public:
  template <typename TCellAccessor>
  inline void
  exchange(
    typename Uni::StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    using CellAccessorType = TCellAccessor;

    using Pointers = Uni::StructuredMemory::Pointers
                     <PetscScalar, CellAccessorType::Dimensions>;

    accessor.pressure() = Pointers::dereference(array, index);
  }
};

class PpeRhsAcquiererAction2 {
public:
  template <typename TCellAccessor>
  inline void
  exchange(
    typename Uni::StructuredMemory::Pointers
    <PetscScalar, TCellAccessor::Dimensions>::Type array,
    typename TCellAccessor::VectorDiType const& index,
    TCellAccessor const& accessor) {
    using CellAccessorType = TCellAccessor;

    using Pointers = Uni::StructuredMemory::Pointers
                     <PetscScalar, CellAccessorType::Dimensions>;

    auto const& value =  Pointers::dereference(array, index);

    accessor.pressure()      += value;
    accessor.projectionTerm() = value;
  }
};
}
}
}
}
