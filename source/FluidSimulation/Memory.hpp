#pragma once

#include <Uni/ExecutionControl/assert>
#include <Uni/ExecutionControl/exception>

#include <array>
#include <memory>
#include <string>

namespace FsiSimulation {
namespace FluidSimulation {
namespace Private {
class Attribute {
public:
  enum class Type {
    Vector = 0,
    Scalar = 1
  };

  Attribute() {}

  Attribute(std::string name_,
            Type const& type_) :
    name(name_),
    type(type_) {}

  Attribute(Attribute const& other) : name(other.name),
    type(other.type) {}

  std::string name;
  Type        type;
};
}

template <typename TSolverTraits>
class Memory {
public:
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

  using AttributeType = Private::Attribute;

  using AttributesType = std::array<Private::Attribute, 2>;

public:
  Memory() {}

  void
  initialize(int const&          rank,
             VectorDiType const& processor_size,
             VectorDiType const& global_cell_size,
             VectorDsType const& geometry_width) {
    _parallelDistribution.initialize(rank,
                                     processor_size,
                                     global_cell_size);
    _gridGeometry.initialize(geometry_width,
                             _parallelDistribution.globalCellSize,
                             _parallelDistribution.corner);

    VectorDiType computational_local_size(
      _parallelDistribution.localCellSize + 2 * VectorDiType::Ones());

    typename GridType::FactoryType cell_accessor_factory =
      [ = ] (BaseGridType const* grid, VectorDiType const& index) {
        return CellAccessorType(this, grid, index);
      };

    _grid.initialize(computational_local_size,
                     cell_accessor_factory);

    _velocity.reset(new VectorDsType[_grid.size().prod()]);
    _pressure.reset(new ScalarType[_grid.size().prod()]);
    _position.reset(new int[_grid.size().prod()]);
    _fgh.reset(new VectorDsType[_grid.size().prod()]);
    _convection.reset(new VectorDsType[_grid.size().prod()]);

    _maxVelocity     = VectorDsType::Zero();
    _dt              = 1.0;
    _time            = 0.0;
    _iterationNumber = 0;

    _attributes[0].name = "velocity";
    _attributes[0].type = AttributeType::Type::Vector;
    _attributes[1].name = "pressure";
    _attributes[1].type = AttributeType::Type::Scalar;

    // logGridInitializationInfo(_grid);
    logParallelTopologyInfo(_parallelDistribution);
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

  ScalarType const&
  attribute(std::size_t const& index,
            short const&       attribute_index,
            unsigned const&    dimension = 0) const {
    return _attribute(index, attribute_index, dimension);
  }

  ScalarType&
  attribute(std::size_t const& index,
            short const&       attribute_index,
            unsigned const&    dimension = 0) {
    return _attribute(index, attribute_index, dimension);
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
  velocity(std::size_t const& index) const {
    return _velocity.get()[index];
  }

  VectorDsType&
  velocity(std::size_t const& index) {
    return _velocity.get()[index];
  }

  ScalarType const&
  velocity(std::size_t const& index, unsigned const& dimension) const {
    return _velocity.get()[index].data()[dimension];
  }

  ScalarType&
  velocity(std::size_t const& index, unsigned const& dimension) {
    return _velocity.get()[index].data()[dimension];
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
  pressure(std::size_t const& index) const {
    return _pressure.get()[index];
  }

  ScalarType&
  pressure(std::size_t const& index) {
    return _pressure.get()[index];
  }

  int const*
  position() const {
    return _position.get();
  }

  int*
  position() {
    return _position.get();
  }

  int const&
  position(std::size_t const& index) const {
    return _position.get()[index];
  }

  int&
  position(std::size_t const& index) {
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
  fgh(std::size_t const& index) const {
    return _fgh.get()[index];
  }

  VectorDsType&
  fgh(std::size_t const& index) {
    return _fgh.get()[index];
  }

  ScalarType const&
  fgh(std::size_t const& index, unsigned const& dimension) const {
    return _fgh.get()[index].data()[dimension];
  }

  ScalarType&
  fgh(std::size_t const& index, unsigned const& dimension) {
    return _fgh.get()[index].data()[dimension];
  }

  VectorDsType*
  convection() {
    return _convection.get();
  }

  VectorDsType const&
  convection(std::size_t const& index) const {
    return _convection.get()[index];
  }

  VectorDsType&
  convection(std::size_t const& index) {
    return _convection.get()[index];
  }

  ScalarType const&
  convection(std::size_t const& index, unsigned const& dimension) const {
    return _convection.get()[index].data()[dimension];
  }

  ScalarType&
  convection(std::size_t const& index, unsigned const& dimension) {
    return _convection.get()[index].data()[dimension];
  }

private:
  ScalarType&
  _attribute(std::size_t const& index,
             short const&       attribute_index,
             int const&         dimension) const {
    switch (attribute_index) {
    case 0:

      return _velocity.get()[index].data()[dimension];

    case 1:

      return _pressure.get()[index];

    default:
      throwException("Invalid attribute index");
    }
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

  std::unique_ptr<VectorDsType> _velocity;
  std::unique_ptr<ScalarType>   _pressure;
  std::unique_ptr<int>          _position;
  std::unique_ptr<VectorDsType> _fgh;
  std::unique_ptr<VectorDsType> _convection;
};
}
}
