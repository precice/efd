#pragma once

#include "FsfdMemory.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
struct FsfdDebugMemoryTraits {
  enum {
    AdditionalAttributeSize = 5
  };
};

template <typename TSolverTraits, typename TDerivedTraits>
class FsfdDebugMemory : public FsfdMemory<TSolverTraits, TDerivedTraits> {
  friend class FsfdMemory<TSolverTraits, TDerivedTraits>;

public:
  using Base = FsfdMemory<TSolverTraits, TDerivedTraits>;

  enum {
    Dimensions    = TSolverTraits::Dimensions,
    AttributeSize = Base::AttributeSize + 5
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
    _diffusion.reset(new VectorDsType[this->_grid.size().prod()]);
    _force.reset(new VectorDsType[this->_grid.size().prod()]);
    _bodyForce.reset(new VectorDsType[this->_grid.size().prod()]);

    this->_attributes[Base::AttributeSize + 0].name
      = "convection";
    this->_attributes[Base::AttributeSize + 0].type
      = Attribute::Type::Vector;
    this->_attributes[Base::AttributeSize + 1].name
      = "diffusion";
    this->_attributes[Base::AttributeSize + 1].type
      = Attribute::Type::Vector;
    this->_attributes[Base::AttributeSize + 2].name
      = "force";
    this->_attributes[Base::AttributeSize + 2].type
      = Attribute::Type::Vector;
    this->_attributes[Base::AttributeSize + 3].name
      = "body-force";
    this->_attributes[Base::AttributeSize + 3].type
      = Attribute::Type::Vector;
    this->_attributes[Base::AttributeSize + 4].name
      = "position-in-respect-to-geometry";
    this->_attributes[Base::AttributeSize + 4].type
      = Attribute::Type::Scalar;
  }

  void
  setDiffusionAt(int const& index, VectorDsType const& value) {
    _diffusion.get()[index] = value;
  }

  void
  setForceAt(int const& index, VectorDsType const& value) {
    _force.get()[index] = value;
  }

  void
  setBodyForceAt(int const& index, VectorDsType const& value) {
    _bodyForce.get()[index] = value;
  }

protected:
  ScalarType
  _attribute(int const& index,
             int const& attribute_index,
             int const& dimension) const {
    switch (attribute_index) {
    case Base::AttributeSize + 0:

      return this->_convection.get()[index].data()[dimension];

    case Base::AttributeSize + 1:

      return this->_diffusion.get()[index].data()[dimension];

    case Base::AttributeSize + 2:

      return _force.get()[index].data()[dimension];

    case Base::AttributeSize + 3:

      return _bodyForce.get()[index].data()[dimension];

    case Base::AttributeSize + 4:

      return this->_position.get()[index];

    default:

      return this->Base::_attribute(index,
                                    attribute_index,
                                    dimension);
    }
  }

  std::unique_ptr<VectorDsType> _diffusion;
  std::unique_ptr<VectorDsType> _force;
  std::unique_ptr<VectorDsType> _bodyForce;
};
}
}
