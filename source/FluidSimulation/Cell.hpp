#ifndef FsiSimulation_FluidSimulation_Cell_hpp
#define FsiSimulation_FluidSimulation_Cell_hpp

#include "BasicCell.hpp"

#include <Eigen/Core>

#include <cassert>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TScalar, int TDimensions>
class Cell;

template <typename TScalar, int TDimensions>
class CellTraits < Cell < TScalar, TDimensions >>
  : public BasicCellTraits < CellTraits < Cell<TScalar, 2 >>> {
  using Base = BasicCellTraits < CellTraits < Cell<TScalar, 2 >>>;

  CellTraits() : Base() {
    logInfo("sdfsdfsdf Gotya2");
  }

public:
  static bool
  initializeAttributes() {
    Base::setAttribute("velocity", 0, Attribute::Type::Vector);
    Base::setAttribute("pressure", 1, Attribute::Type::Scalar);

    return true;
  }
};

template <typename TScalar, int TDimensions>
class Cell {
  static_assert((TDimensions > 1) && (TDimensions < 4),
                "Only 2D and 3D simulations supported");

public:
  typedef TScalar                               Scalar;
  typedef Eigen::Matrix<Scalar, TDimensions, 1> VectorDs;
  typedef Eigen::Matrix<int, TDimensions, 1>    VectorDi;
  typedef VectorDs                              Velocity;
  typedef Scalar                                Pressure;

  enum {
    Dimensions = TDimensions
  };

  using Traits = CellTraits < Cell < Scalar, Dimensions >>;

public:
  Cell() {}

  Cell(Cell const& other) = delete;

  ~Cell() {}

  Cell&
  operator=(Cell const& other) = delete;

  Scalar&
  attribute(int const& index, int dimension = 0) {
    switch (index) {
    case 0:

      return _velocity(dimension);

    case 1:

      return _pressure;
    }

    return _pressure;
  }
  Scalar const&
  attribute(int const& index, int dimension = 0) const {
    switch (index) {
    case 0:

      return _velocity(dimension);

    case 1:

      return _pressure;
    }

    return _pressure;
  }

  VectorDs&
  velocity() {
    return _velocity;
  }

  VectorDs const&
  velocity() const {
    return _velocity;
  }

  Scalar&
  velocity(int const& i) {
    return _velocity.data()[i];
  }

  Scalar const&
  velocity(int const& i) const {
    return _velocity.data()[i];
  }

  VectorDs&
  fgh() {
    return _fgh;
  }

  VectorDs const&
  fgh() const {
    return _fgh;
  }

  Scalar&
  fgh(int const& i) {
    return _fgh.data()[i];
  }

  Scalar const&
  fgh(int const& i) const {
    return _fgh.data()[i];
  }

  VectorDs&
  convection() {
    return _convection;
  }

  VectorDs const&
  convection() const {
    return _convection;
  }

  Scalar&
  convection(int const& i) {
    return _convection.data()[i];
  }

  Scalar const&
  convection(int const& i) const {
    return _convection.data()[i];
  }

  Scalar&
  pressure() {
    return _pressure;
  }

  Scalar const&
  pressure() const {
    return _pressure;
  }

  Scalar&
  pressureProjection() {
    return _pressureProjection;
  }

  Scalar const&
  pressureProjection() const {
    return _pressureProjection;
  }

  int&
  position() {
    return _position;
  }

  int const&
  position() const {
    return _position;
  }

  VectorDi&
  positions() {
    return _positions;
  }

  int&
  positions(int const& dimension) {
    return _positions(dimension);
  }

  VectorDi const&
  positions() const {
    return _positions;
  }

  int const&
  positions(int const& dimension) const {
    return _positions(dimension);
  }

private:
  VectorDs _velocity;
  VectorDs _fgh;
  VectorDs _convection;
  Scalar   _pressure;
  Scalar   _pressureProjection;
  int      _position;
  VectorDi _positions;

  static Traits const traits;
};

template <typename TScalar, int TDimensions>
typename Cell<TScalar, TDimensions>::Traits const Cell<TScalar,
                                                       TDimensions>::traits;
}
}

#endif
