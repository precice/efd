#pragma once

#include "SfsfdCellAccessor.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class IfsfdCellAccessor : public SfsfdCellAccessor<TSolverTraits> {
public:
  using BaseType = SfsfdCellAccessor<TSolverTraits>;

  enum {
    Dimensions = TSolverTraits::Dimensions
  };
  using MemoryType = typename TSolverTraits::MemoryType;

  using GridType = typename TSolverTraits::GridType;

  using BaseGridType = typename TSolverTraits::BaseGridType;

  using VectorDsType = typename TSolverTraits::VectorDsType;

  using VectorDiType = typename TSolverTraits::VectorDiType;

  using ScalarType = typename TSolverTraits::ScalarType;

public:
  IfsfdCellAccessor(MemoryType*         memory,
                    BaseGridType const* grid) : BaseType(memory, grid) {}

  IfsfdCellAccessor(MemoryType*         memory,
                    BaseGridType const* grid,
                    VectorDiType const& index) : BaseType(memory, grid,
                                                          index) {}

  IfsfdCellAccessor(IfsfdCellAccessor const& other) : BaseType(other) {}

  ~IfsfdCellAccessor() {}

  IfsfdCellAccessor&
  operator=(IfsfdCellAccessor const& other) {
    this->BaseType::operator=(other);

    return *this;
  }

  ScalarType&
  projectionTerm() const {
    return this->_memory->projectionTerm(this->globalIndex());
  }

  ScalarType&
  projectionTerm(int const& dimension, int const& offset) const {
    return this->_memory->projectionTerm(
      this->relativeGlobalIndex(dimension, offset));
  }

  ScalarType&
  projectionTerm(int const& dimension,
                 int const& offset,
                 int const& dimension2,
                 int const& offset2) const {
    return this->_memory->projectionTerm(
      this->relativeGlobalIndex(dimension, offset,
                                dimension2, offset2));
  }

  ScalarType&
  projectionTerm(VectorDiType const& index) const {
    return this->_memory->projectionTerm(this->relativeGlobalIndex(index));
  }
};
}
}
