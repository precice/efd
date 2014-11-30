#ifndef FsiSimulation_Cell_hpp
#define FsiSimulation_Cell_hpp

#include <Eigen/Core>

namespace FsiSimulation {
template <typename Scalar, int D>
class Cell {
public:
  typedef Eigen::Matrix<Scalar, D, 1> VectorDs;

public:
  Cell() {}

  Cell(Cell const& other) = delete;

  ~Cell() {}

  Cell&
  operator=(Cell const& other) = delete;

  VectorDs&
  velocity() { return _velocity; }

  VectorDs const&
  velocity() const { return _velocity; }

  Scalar&
  velocity(int const& i) { return _velocity.data()[i]; }

  Scalar const&
  velocity(int const& i) const { return _velocity.data()[i]; }

  VectorDs&
  fgh() { return _fgh; }

  VectorDs const&
  fgh() const { return _fgh; }

  Scalar&
  fgh(int const& i) { return _fgh.data()[i]; }

  Scalar const&
  fgh(int const& i) const { return _fgh.data()[i]; }

  Scalar&
  pressure() { return _pressure; }

  Scalar const&
  pressure() const { return _pressure; }

private:
  VectorDs _velocity;
  VectorDs _fgh;
  Scalar   _pressure;
};
}

#endif
