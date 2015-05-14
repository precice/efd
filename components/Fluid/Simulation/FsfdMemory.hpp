#pragma once

#include "Grid.hpp"
#include "GridGeometry.hpp"
#include "ParallelDistribution.hpp"
#include "Parameters.hpp"

#include <Uni/ExecutionControl/assert>
#include <Uni/ExecutionControl/exception>

#include <array>
#include <memory>
#include <string>

namespace FsiSimulation {
namespace FluidSimulation {
class Attribute {
public:
  enum class Type {
    Vector,
    Scalar
  };

  enum class DisplayMode {
    Cells,
    Points
  };

  Attribute() :
    mode(DisplayMode::Cells),
    doCentralize(false) {}

  Attribute(std::string        name_,
            Type const&        type_,
            DisplayMode const& mode_)
    : name(name_),
    type(type_),
    mode(mode_),
    doCentralize(false) {}

  Attribute(Attribute const& other)
    : name(other.name),
    type(other.type),
    mode(other.mode),
    doCentralize(other.doCentralize) {}

  std::string name;
  Type        type;
  DisplayMode mode;
  bool        doCentralize;
};

namespace Private {
template <unsigned TDebugLeve>
struct AttributesSizeTraits {};

template <>
struct AttributesSizeTraits<0> {
  enum {
    Size = 2
  };
};

template <>
struct AttributesSizeTraits<1> {
  enum {
    Size = AttributesSizeTraits<0>::Size + 5
  };
};
}

template <typename TSolverTraits>
class FsfdMemory {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions     = TSolverTraits::Dimensions,
    DebugLevel     = SolverTraitsType::DebugLevel,
    AttributesSize = Private::AttributesSizeTraits<DebugLevel>::Size
  };

  using GridGeometryType = typename TSolverTraits::GridGeometryType;

  using ParallelDistributionType
          = typename TSolverTraits::ParallelDistributionType;

  using ParametersType = typename TSolverTraits::ParametersType;

  using GridType = typename TSolverTraits::GridType;

  using SubgridType = typename TSolverTraits::SubgridType;

  using CellAccessorType = typename TSolverTraits::CellAccessorType;

  using MemoryType = typename TSolverTraits::MemoryType;

  using VectorDsType = typename TSolverTraits::VectorDsType;

  using VectorDiType = typename TSolverTraits::VectorDiType;

  using ScalarType = typename TSolverTraits::ScalarType;

  using LocationType = Eigen::Matrix<int, Dimensions + 1, 1>;

  using AttributeType = Attribute;

  using AttributesType = std::array<Attribute, AttributesSize>;

public:
  FsfdMemory() {}

  void
  initialize(VectorDiType const& processor_size,
             VectorDiType const& global_cell_size,
             VectorDsType const& geometry_width) {
    _parallelDistribution.initialize(processor_size,
                                     global_cell_size,
                                     VectorDiType::Ones());
    _gridGeometry.initialize(geometry_width,
                             _parallelDistribution.globalCellSize,
                             _parallelDistribution.corner);

    VectorDiType computational_local_size(
      _parallelDistribution.localCellSize + 2 * VectorDiType::Ones());

    typename GridType::FactoryType cell_accessor_factory
      = [ = ] (SubgridType const* grid, VectorDiType const& index) {
          return CellAccessorType(static_cast<MemoryType*>(this),
                                  grid,
                                  index);
        };

    _grid.initialize(computational_local_size,
                     cell_accessor_factory);

    _velocity.reset(new VectorDsType[_grid.size().prod()]);
    _pressure.reset(new ScalarType[_grid.size().prod()]);
    _position.reset(new LocationType[_grid.size().prod()]);
    _fgh.reset(new VectorDsType[_grid.size().prod()]);
    _convection.reset(new VectorDsType[_grid.size().prod()]);

    _maxVelocity     = VectorDsType::Zero();
    _dt              = 1.0;
    _time            = 0.0;
    _iterationNumber = 0;

    _attributes[0].name         = "velocity";
    _attributes[0].type         = AttributeType::Type::Vector;
    // _attributes[0].mode         = AttributeType::DisplayMode::Points;
    _attributes[0].doCentralize = true;
    _attributes[1].name         = "pressure";
    _attributes[1].type         = AttributeType::Type::Scalar;
    // _attributes[1].mode         = AttributeType::DisplayMode::Points;

    if (DebugLevel == 1) {
      _diffusion.reset(new VectorDsType[this->_grid.size().prod()]);
      _force.reset(new VectorDsType[this->_grid.size().prod()]);
      _bodyForce.reset(new VectorDsType[this->_grid.size().prod()]);
      this->_attributes[2].name
        = "convection";
      this->_attributes[2].type
        = Attribute::Type::Vector;
      this->_attributes[3].name
        = "diffusion";
      this->_attributes[3].type
        = Attribute::Type::Vector;
      this->_attributes[4].name
        = "force";
      this->_attributes[4].type
        = Attribute::Type::Vector;
      this->_attributes[5].name
        = "body-force";
      this->_attributes[5].type
        = Attribute::Type::Vector;
      this->_attributes[6].name
        = "position-in-respect-to-geometry";
      this->_attributes[6].type
        = Attribute::Type::Scalar;
    }

    // logGridInitializationInfo(_grid);
    _parallelDistribution.toString();
    logInfo(_gridGeometry.toString());
  }

  ParallelDistributionType const*
  parallelDistribution() const {
    return &_parallelDistribution;
  }

  ParallelDistributionType*
  parallelDistribution() {
    return &_parallelDistribution;
  }

  GridGeometryType const*
  gridGeometry() const {
    return &_gridGeometry;
  }

  GridGeometryType*
  gridGeometry() {
    return &_gridGeometry;
  }

  GridType const*
  grid() const {
    return &_grid;
  }

  GridType*
  grid() {
    return &_grid;
  }

  ParametersType const*
  parameters() const {
    return &_parameters;
  }

  ParametersType*
  parameters() {
    return &_parameters;
  }

  AttributesType const*
  attributes() const {
    return &_attributes;
  }

  AttributesType*
  attributes() {
    return &_attributes;
  }

  VectorDsType const&
  maxVelocity() const {
    return _maxVelocity;
  }

  VectorDsType&
  maxVelocity() {
    return _maxVelocity;
  }

  ScalarType const&
  timeStepSize() const {
    return _dt;
  }

  ScalarType&
  timeStepSize() {
    return _dt;
  }

  long double const&
  time() const {
    return _time;
  }

  long double&
  time() {
    return _time;
  }

  unsigned long long const&
  iterationNumber() const {
    return _iterationNumber;
  }

  unsigned long long&
  iterationNumber() {
    return _iterationNumber;
  }

  // Those members are for the debug reasons {{{
  void
  setDiffusionAt(int const& index, VectorDsType const& value) {
    if (DebugLevel > 0) {
      _diffusion.get()[index] = value;
    }
  }

  void
  setForceAt(int const& index, VectorDsType const& value) {
    if (DebugLevel > 0) {
      _force.get()[index] = value;
    }
  }

  void
  addForceAt(int const& index, VectorDsType const& value) {
    if (DebugLevel > 0) {
      _force.get()[index] += value;
    }
  }

  void
  setBodyForceAt(int const& index, VectorDsType const& value) {
    if (DebugLevel > 0) {
      _bodyForce.get()[index] = value;
    }
  }
  // Those members are for the debug reasons }}}

  ScalarType
  attribute(int const&   index,
            short const& attribute_index,
            int const&   dimension = 0) const {
    return static_cast<MemoryType const*>(this)
           ->_attribute(index,
                        attribute_index,
                        dimension);
  }

  VectorDsType const*
  velocity() const {
    return _velocity.get();
  }

  VectorDsType*
  velocity() {
    return _velocity.get();
  }

  VectorDsType const&
  velocity(int const& index) const {
    return _velocity.get()[index];
  }

  VectorDsType&
  velocity(int const& index) {
    return _velocity.get()[index];
  }

  ScalarType const&
  velocity(int const& index, int const& dimension) const {
    return _velocity.get()[index](dimension);
  }

  ScalarType&
  velocity(int const& index, int const& dimension) {
    return _velocity.get()[index](dimension);
  }

  ScalarType const*
  pressure() const {
    return _pressure.get();
  }

  ScalarType*
  pressure() {
    return _pressure.get();
  }

  ScalarType const&
  pressure(int const& index) const {
    return _pressure.get()[index];
  }

  ScalarType&
  pressure(int const& index) {
    return _pressure.get()[index];
  }

  LocationType const*
  position() const {
    return _position.get();
  }

  LocationType*
  position() {
    return _position.get();
  }

  LocationType const&
  position(int const& index) const {
    return _position.get()[index];
  }

  LocationType&
  position(int const& index) {
    return _position.get()[index];
  }

  VectorDsType const*
  fgh() const {
    return _fgh.get();
  }

  VectorDsType*
  fgh() {
    return _fgh.get();
  }

  VectorDsType const&
  fgh(int const& index) const {
    return _fgh.get()[index];
  }

  VectorDsType&
  fgh(int const& index) {
    return _fgh.get()[index];
  }

  ScalarType const&
  fgh(int const& index, int const& dimension) const {
    return _fgh.get()[index].data()[dimension];
  }

  ScalarType&
  fgh(int const& index, int const& dimension) {
    return _fgh.get()[index].data()[dimension];
  }

  VectorDsType*
  convection() {
    return _convection.get();
  }

  VectorDsType const&
  convection(int const& index) const {
    return _convection.get()[index];
  }

  VectorDsType&
  convection(int const& index) {
    return _convection.get()[index];
  }

  ScalarType const&
  convection(int const& index, int const& dimension) const {
    return _convection.get()[index].data()[dimension];
  }

  ScalarType&
  convection(int const& index, int const& dimension) {
    return _convection.get()[index].data()[dimension];
  }

protected:
  ScalarType
  _attribute(int const& index,
             int const& attribute_index,
             int const& dimension) const {
    if (DebugLevel > 0) {
      switch (attribute_index) {
      case 2:

        return this->_convection.get()[index](dimension);

      case 3:

        return this->_diffusion.get()[index](dimension);

      case 4:

        return _force.get()[index](dimension);

      case 5:

        return _bodyForce.get()[index](dimension);

      case 6:

        return static_cast<ScalarType>(
          this->_position.get()[index](Dimensions));
      }
    }

    switch (attribute_index) {
    case 0:

      return _velocity.get()[index].data()[dimension];

    case 1:

      return _pressure.get()[index];

    default:
      throwException("Invalid attribute index");
    }

    return 0.0;
  }

  ParallelDistributionType _parallelDistribution;
  GridGeometryType         _gridGeometry;
  GridType                 _grid;
  ParametersType           _parameters;
  AttributesType           _attributes;
  VectorDsType             _maxVelocity;
  ScalarType               _dt;
  long double              _time;
  unsigned long long       _iterationNumber;

  std::unique_ptr<VectorDsType[]> _velocity;
  std::unique_ptr<ScalarType[]>   _pressure;
  std::unique_ptr<LocationType>   _position;
  std::unique_ptr<VectorDsType[]> _fgh;
  std::unique_ptr<VectorDsType[]> _convection;

  // The members below are for the debug reasons
  std::unique_ptr<VectorDsType[]> _diffusion;
  std::unique_ptr<VectorDsType[]> _force;
  std::unique_ptr<VectorDsType[]> _bodyForce;
};
}
}
