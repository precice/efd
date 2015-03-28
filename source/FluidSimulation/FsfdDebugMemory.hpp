#pragma once

#include "FsfdMemory.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
struct FsfdDebugMemoryTraits {
  enum {
    AdditionalAttributeSize = 3
  };
};

template <typename TSolverTraits, typename TDerivedTraits>
class FsfdDebugMemory : public FsfdMemory<TSolverTraits, TDerivedTraits> {
  friend class FsfdMemory<TSolverTraits, TDerivedTraits>;

public:
  using Base = FsfdMemory<TSolverTraits, TDerivedTraits>;

  enum {
    Dimensions    = TSolverTraits::Dimensions,
    AttributeSize = Base::AttributeSize + 3
  };

  using VectorDsType = typename TSolverTraits::VectorDsType;

  using VectorDiType = typename TSolverTraits::VectorDiType;

  using ScalarType = typename TSolverTraits::ScalarType;

public:
  FsfdDebugMemory() {}

  void
  initialize(VectorDiType const& processor_size,
             VectorDiType const& global_cell_size,
             VectorDsType const& geometry_width) {
    this->Base::initialize(processor_size,
                           global_cell_size,
                           geometry_width);
    _force.reset(new VectorDsType[this->_grid.size().prod()]);
    _bodyForce.reset(new VectorDsType[this->_grid.size().prod()]);
    this->_attributes[Base::AttributeSize + 0].name
      = "force";
    this->_attributes[Base::AttributeSize + 0].type
      = Attribute::Type::Vector;
    this->_attributes[Base::AttributeSize + 1].name
      = "body-force";
    this->_attributes[Base::AttributeSize + 1].type
      = Attribute::Type::Vector;
    this->_attributes[Base::AttributeSize + 2].name
      = "position-in-respect-to-geometry";
    this->_attributes[Base::AttributeSize + 2].type
      = Attribute::Type::Scalar;
  }

  void
  setForceAt(int const& index, VectorDsType const& force) {
    _force.get()[index] = force;
  }

  void
  setBodyForceAt(int const& index, VectorDsType const& force) {
    _bodyForce.get()[index] = force;
  }

protected:
  ScalarType
  _attribute(int const& index,
             int const& attribute_index,
             int const& dimension) const {
    switch (attribute_index) {
    case Base::AttributeSize + 0:

      return _force.get()[index].data()[dimension];

    case Base::AttributeSize + 1:

      return _bodyForce.get()[index].data()[dimension];

    case Base::AttributeSize + 2:

      return this->_position.get()[index];

    default:

      return this->Base::_attribute(index,
                                    attribute_index,
                                    dimension);
    }
  }

  std::unique_ptr<VectorDsType> _force;
  std::unique_ptr<VectorDsType> _bodyForce;
};
}
}
