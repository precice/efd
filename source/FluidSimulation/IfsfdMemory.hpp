#pragma once

#include "SfsfdMemory.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class IfsfdMemory : public SfsfdMemory<TSolverTraits> {
public:
  using BaseType = SfsfdMemory<TSolverTraits>;

  enum {
    Dimensions = TSolverTraits::Dimensions
  };

  using GridGeometryType = typename TSolverTraits::GridGeometryType;

  using ParallelDistributionType
          = typename TSolverTraits::ParallelDistributionType;

  using ParametersType = typename TSolverTraits::ParametersType;

  using GridType = typename TSolverTraits::GridType;

  using BaseGridType = typename TSolverTraits::BaseGridType;

  using CellAccessorType = typename TSolverTraits::CellAccessorType;

  using VectorDsType = typename TSolverTraits::VectorDsType;

  using VectorDiType = typename TSolverTraits::VectorDiType;

  using ScalarType = typename TSolverTraits::ScalarType;

public:
  IfsfdMemory() {}

  void
  initialize(VectorDiType const& processor_size,
             VectorDiType const& global_cell_size,
             VectorDsType const& geometry_width) {
    this->BaseType::initialize(processor_size,
                               global_cell_size,
                               geometry_width);

    _projectionTerm.reset(new ScalarType[this->grid()->size().prod()]);
  }

  ScalarType const*
  projectionTerm() const {
    return _projectionTerm.get();
  }

  ScalarType*
  projectionTerm() {
    return _projectionTerm.get();
  }

  ScalarType const&
  projectionTerm(std::size_t const& index) const {
    return _projectionTerm.get()[index];
  }

  ScalarType&
  projectionTerm(std::size_t const& index) {
    return _projectionTerm.get()[index];
  }

private:
  std::unique_ptr<ScalarType> _projectionTerm;
};
}
}
