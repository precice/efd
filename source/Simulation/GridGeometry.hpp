#pragma once

#include <Uni/Logging/format>

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TScalar, unsigned TDimensions>
class UniformGridGeometry {
public:
  using Scalar = TScalar;

  enum {
    Dimensions = TDimensions
  };

  using VectorDi =  Eigen::Matrix<int, Dimensions, 1>;

  using VectorDs = Eigen::Matrix<Scalar, Dimensions, 1>;

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

  Scalar const&
  size(int const& index) const {
    return _size(index);
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

  VectorDs
  computeCellPosition(VectorDi const& i) const {
    return _cellWidth.cwiseProduct(i.template cast<Scalar>());
  }

  Scalar
  computeCellPosition(int const& dimension, int const& i) const {
    return _cellWidth(dimension) * i;
  }

  VectorDs const&
  minCellWidth() const {
    return _cellWidth;
  }

  Scalar const&
  minCellWidth(unsigned const& dimension) const {
    return _cellWidth(dimension);
  }

  VectorDs const&
  maxCellWidth() const {
    return _cellWidth;
  }

  Scalar&
  minCellWidth(unsigned const& dimension) {
    return _cellWidth(dimension);
  }

  std::string
  toString() {
    return Uni::Logging::format("{1}\n"
                                "{2}\n"
                                "{3}\n",
                                _size.transpose(),
                                _cellWidth.transpose(),
                                _corner.transpose());
  }

private:
  VectorDs _size;
  VectorDs _cellWidth;
  VectorDi _corner;
};
}
}
