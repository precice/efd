#ifndef Fluid_Simulation_IfsfdMemory_hpp
#define Fluid_Simulation_IfsfdMemory_hpp

#include "FsfdMemory.hpp"

namespace Fluid {
namespace Simulation {
template <typename TSolverTraits>
class IfsfdMemory : public FsfdMemory<TSolverTraits> {
  friend class FsfdMemory<TSolverTraits>;

public:
  using BaseType = FsfdMemory<TSolverTraits>;

  enum {
    Dimensions     = TSolverTraits::Dimensions,
    AttributesSize = BaseType::AttributesSize
  };

  using GridGeometryType = typename TSolverTraits::GridGeometryType;

  using ParallelDistributionType
          = typename TSolverTraits::ParallelDistributionType;

  using ParametersType = typename TSolverTraits::ParametersType;

  using GridType = typename TSolverTraits::GridType;

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

protected:
  ScalarType
  _attribute(int const& index,
             int const& attribute_index,
             int const& dimension) const {
    return this->BaseType::_attribute(index, attribute_index, dimension);
  }

private:
  std::unique_ptr<ScalarType[]> _projectionTerm;
};
}
}
#endif
