#ifndef FsiSimulation_FluidSimulation_GridGeomtry_hpp
#define FsiSimulation_FluidSimulation_GridGeomtry_hpp

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TScalar, int TD>
class UniformGridGeometry {
public:
  typedef Eigen::Matrix<int, TD, 1>     VectorDi;
  typedef Eigen::Matrix<TScalar, TD, 1> VectorDs;

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
    _cellWidth = _size.cwiseQuotient(cellSize.template cast<TScalar>());
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
    return _cellWidth.cwiseProduct((_corner + i).template cast<TScalar>());
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
}
#endif
