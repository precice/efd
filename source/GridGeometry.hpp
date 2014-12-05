#ifndef FsiSimulation_GridGeomtry
#define FsiSimulation_GridGeomtry

#include <Eigen/Core>

namespace FsiSimulation {
template <typename Scalar, int D>
class UniformGridGeometry {
public:
  typedef Eigen::Matrix<int, D, 1>    VectorDi;
  typedef Eigen::Matrix<Scalar, D, 1> VectorDs;

public:
  UniformGridGeometry() {}

  UniformGridGeometry(UniformGridGeometry const& other) = delete;

  ~UniformGridGeometry() {}

  UniformGridGeometry&
  operator=(UniformGridGeometry const& other) = delete;

  void
  initialize(VectorDs const& size,
             VectorDi const& cellSize,
             VectorDi const& corner) {
    _size      = size;
    _cellWidth = _size.cwiseQuotient(cellSize.template cast<Scalar>());
    _corner    = corner;
  }

  VectorDs const&
  size() const {
    return _size;
  }

  VectorDs const&
  cellWidth(VectorDi const& i) const {
    ((void)i);

    return _cellWidth;
  }

  VectorDi const&
  corner() const {
    return _corner;
  }

  VectorDs
  cellPosition(VectorDi const& i) const {
    return _cellWidth.cwiseProduct((_corner + i).template cast<Scalar>());
  }

  VectorDs const&
  minCellWidth() const {
    return _cellWidth;
  }

private:
  VectorDs _size;
  VectorDs _cellWidth;
  VectorDi _corner;
};
}
#endif
