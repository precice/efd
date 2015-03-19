#pragma once

#include "BasicCellAccessor.hpp"
#include "Grid.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class CellAccessor;
template <typename TSolverTraits>
class CellAccessor : public BasicCellAccessor<TSolverTraits> {
public:
  using BaseType = BasicCellAccessor<TSolverTraits>;

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
  CellAccessor(MemoryType *   memory,
               BaseGridType const* grid) :
    BaseType(memory, grid) {}

  CellAccessor(MemoryType *   memory,
               BaseGridType const* grid,
               VectorDiType const& index)
    : BaseType(memory, grid, index) {}

  CellAccessor(CellAccessor const& other)
    : BaseType(other) {}

  ~CellAccessor() {}

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->BaseType::operator=(other);

    return *this;
  }
};
}
}
