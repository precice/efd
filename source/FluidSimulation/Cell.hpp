#ifndef FsiSimulation_FluidSimulation_Cell_hpp
#define FsiSimulation_FluidSimulation_Cell_hpp

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TScalar, int TD>
class Cell {
public:
  typedef TScalar                      Scalar;
  typedef Eigen::Matrix<Scalar, TD, 1> VectorDs;
  typedef Eigen::Matrix<int, TD, 1>    VectorDi;
  typedef VectorDs                     Velocity;
  typedef Scalar                       Pressure;

  enum {
    Dimensions = TD
  };

public:
  Cell() {}

  Cell(Cell const& other) = delete;

  ~Cell() {}

  Cell&
  operator=(Cell const& other) = delete;

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

  VectorDs&
  distances() {
    return _distances;
  }

  VectorDs const&
  distances() const {
    return _distances;
  }

  Scalar&
  distances(int const& i) {
    return _distances.data()[i];
  }

  Scalar const&
  distances(int const& i) const {
    return _distances.data()[i];
  }

  Scalar&
  pressure() {
    return _pressure;
  }

  Scalar const&
  pressure() const {
    return _pressure;
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
  VectorDs _distances;
  Scalar   _pressure;
  int _position;
  VectorDi _positions;
};
}
}

#endif
