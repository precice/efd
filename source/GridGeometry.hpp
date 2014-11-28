#ifndef FsiSimulation_GridGeomtry
#define FsiSimulation_GridGeomtry

#include <Eigen/Core>

namespace FsiSimulation {
template <typename Scalar, int D>
class UniformGridGeometry {
public:
  typedef Eigen::Matrix<Scalar, D, 1> VectorDi;
  typedef Eigen::Matrix<Scalar, D, 1> VectorDs;

public:
  UniformGridGeometry() {}

  UniformGridGeometry(UniformGridGeometry const& other) = delete;

  ~UniformGridGeometry() {}

  UniformGridGeometry&
  operator=(UniformGridGeometry const& other) = delete;

  VectorDs
  cellSize(VectorDi const& i) const {
    ((void)i);

    return _cellSize;
  }

  void
  cellSize(VectorDs const& cellSize) {
    _cellSize = cellSize;
  }

  void
  corner(VectorDi const& corner) {
    _corner = corner;
  }

  VectorDs
  cellPosition(VectorDi const& i) const {
    return _cellSize * (_corner + i);
  }

  VectorDs
  cellMinSize() const {
    return _cellSize;
  }

private:
  VectorDs _cellSize;
  VectorDi _corner;
};
}
#endif
