#pragma once

#include <Uni/StructuredGrid/Basic/GlobalMultiIndex>

namespace FsiSimulation {
namespace FluidSimulation {
namespace Private {
template <typename TSolverTraits>
struct MultiIndexTraits {
  using Type = typename TSolverTraits::CellAccessorType;

  using GridType = typename TSolverTraits::BaseGridType;

  enum {
    Dimensions = TSolverTraits::Dimensions
  };
};
}
template <typename TSolverTraits>
class BasicCellAccessor
  : public Uni::StructuredGrid::Basic::GlobalMultiIndex
    < Private::MultiIndexTraits < TSolverTraits >> {
public:
  using BaseType = Uni::StructuredGrid::Basic::GlobalMultiIndex
                   < Private::MultiIndexTraits < TSolverTraits >>;

  enum {
    Dimensions = TSolverTraits::Dimensions
  };

  using MemoryType = typename TSolverTraits::MemoryType;

  using GridType = typename TSolverTraits::GridType;

  using BaseGridType = typename TSolverTraits::BaseGridType;

  using VectorDsType = typename TSolverTraits::VectorDsType;

  using VectorDiType = typename TSolverTraits::VectorDiType;

  using ScalarType = typename TSolverTraits::ScalarType;

public:
  BasicCellAccessor(MemoryType*         memory,
                    BaseGridType const* grid)
    : BaseType(grid),
    _memory(memory) {}

  BasicCellAccessor(MemoryType*         memory,
                    BaseGridType const* grid,
                    VectorDiType const& index)
    : BaseType(grid, index),
    _memory(memory) {}

  BasicCellAccessor(BasicCellAccessor const& other)
    : BaseType(other),
    _memory(other._memory) {}

  ~BasicCellAccessor() {}

  BasicCellAccessor&
  operator=(BasicCellAccessor const& other) {
    this->BaseType::operator=(other);
    _memory = other._memory;

    return *this;
  }

  MemoryType*
  memory() const {
    return _memory;
  }

  ScalarType&
  attribute(int const& attribute_index,
            int const& dimension = 0) const {
    return _attribute(this->globalIndex(), attribute_index, dimension);
  }

  VectorDsType
  width() const {
    return _memory->gridGeometry()->cellWidth(
      this->index() - _memory->grid()->leftIndent());
  }

  VectorDsType
  width(int const& dimension, int const& offset) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(dimension, offset)
      - _memory->grid()->leftIndent());
  }

  VectorDsType
  width(int const& dimension,
        int const& offset,
        int const& dimension2,
        int const& offset2) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(dimension,
                          offset,
                          dimension2,
                          offset2)
      - _memory->grid()->leftIndent());
  }

  VectorDsType
  width(VectorDiType const& index) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(index) - _memory->grid()->leftIndent());
  }

  ScalarType
  width(int const& dimension) const {
    return _memory->gridGeometry()->cellWidth(
      this->indexValues() - _memory->grid()->leftIndent())(dimension);
  }

  ScalarType
  width(int const& dimension,
        int const& offset,
        int const& width_dimension) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(dimension, offset) - _memory->grid()->leftIndent())
             (width_dimension);
  }

  ScalarType
  width(int const& dimension,
        int const& offset,
        int const& dimension2,
        int const& offset2,
        int const& width_dimension) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(dimension,
                          offset,
                          dimension2,
                          offset2)
      - _memory->grid()->leftIndent())
             (width_dimension);
  }

  VectorDsType
  width(VectorDiType const& index,
        int const&          width_dimension) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(index) - _memory->grid()->leftIndent())
             (width_dimension);
  }

  VectorDsType
  position() const {
    return _memory->gridGeometry()->cellPosition(
      this->index() - _memory->grid()->leftIndent());
  }

  VectorDsType
  position(int const& dimension, int const& offset) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(dimension, offset) - _memory->grid()->leftIndent());
  }

  VectorDsType
  position(int const& dimension,
           int const& offset,
           int const& dimension2,
           int const& offset2) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(dimension,
                          offset,
                          dimension2,
                          offset2)
      - _memory->grid()->leftIndent());
  }

  VectorDsType
  position(VectorDiType const& index) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(index)
      - _memory->grid()->leftIndent());
  }

  ScalarType
  position(int const& dimension) const {
    return _memory->gridGeometry()->cellPosition(
      this->indexValues() - _memory->grid()->leftIndent())(dimension);
  }

  ScalarType
  position(int const& dimension,
           int const& offset,
           int const& position_dimension) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(dimension, offset) - _memory->grid()->leftIndent())
             (position_dimension);
  }

  ScalarType
  position(int const& dimension,
           int const& offset,
           int const& dimension2,
           int const& offset2,
           int const& position_dimension) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(dimension,
                          offset,
                          dimension2,
                          offset2)
      - _memory->grid()->leftIndent())
             (position_dimension);
  }

  VectorDsType
  position(VectorDiType const& index,
           int const&          position_dimension) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(index) - _memory->grid()->leftIndent());
  }

  VectorDsType
  velocityPosition(int const& velocity_dimension) const {
    return _velocityPosition(this->index(), velocity_dimension);
  }

  VectorDsType
  velocityPosition(int const& dimension,
                   int const& offset,
                   int const& velocity_dimension) const {
    VectorDiType index = this->relativeIndex(dimension, offset);

    return _velocityPosition(index, velocity_dimension);
  }

  VectorDsType
  velocityPosition(int const& dimension,
                   int const& offset,
                   int const& dimension2,
                   int const& offset2,
                   int const& velocity_dimension) const {
    VectorDiType index = this->relativeIndex(dimension, offset,
                                             dimension2, offset2);

    return _velocityPosition(index, velocity_dimension);
  }

  VectorDsType
  velocityPosition(VectorDiType const& relative_index,
                   int const&          velocity_dimension) const {
    VectorDiType index = this->relativeIndex(relative_index);

    return _velocityPosition(index, velocity_dimension);
  }

  ScalarType
  velocityPosition(int const& velocity_dimension,
                   int const& position_dimension) const {
    return _velocityPosition(this->index(),
                             velocity_dimension,
                             position_dimension);
  }

  ScalarType
  velocityPosition(int const& dimension,
                   int const& offset,
                   int const& velocity_dimension,
                   int const& position_dimension) const {
    VectorDiType index = this->relativeIndex(dimension, offset);

    return _velocityPosition(index,
                             velocity_dimension,
                             position_dimension);
  }

  ScalarType
  velocityPosition(int const& dimension,
                   int const& offset,
                   int const& dimension2,
                   int const& offset2,
                   int const& velocity_dimension,
                   int const& position_dimension) const {
    VectorDiType index = this->relativeIndex(dimension, offset,
                                             dimension2, offset2);

    return _velocityPosition(index,
                             velocity_dimension,
                             position_dimension);
  }

  ScalarType
  velocityPosition(VectorDiType const& relative_index,
                   int const&          velocity_dimension,
                   int const&          position_dimension) const {
    VectorDiType index = this->relativeIndex(relative_index);

    return _velocityPosition(index,
                             velocity_dimension,
                             position_dimension);
  }

  VectorDsType
  pressurePosition() const {
    return _pressurePosition(this->index());
  }

  VectorDsType
  pressurePosition(int const& dimension,
                   int const& offset) const {
    return _pressurePosition(this->relativeIndex(dimension, offset));
  }

  VectorDsType
  pressurePosition(int const& dimension,
                   int const& offset,
                   int const& dimension2,
                   int const& offset2) const {
    VectorDiType index = this->relativeIndex(dimension, offset,
                                             dimension2, offset2);

    return _pressurePosition(index);
  }

  VectorDsType
  pressurePosition(VectorDiType const& relative_index) const {
    VectorDiType index = this->relativeIndex(relative_index);

    return _pressurePosition(index);
  }

  ScalarType
  pressurePosition(int const& position_dimension) const {
    return _pressurePosition(this->index(), position_dimension);
  }

  ScalarType
  pressurePosition(int const& dimension,
                   int const& offset,
                   int const& position_dimension) const {
    return _pressurePosition(this->relativeIndex(dimension, offset),
                             position_dimension);
  }

  ScalarType
  pressurePosition(int const& dimension,
                   int const& offset,
                   int const& dimension2,
                   int const& offset2,
                   int const& position_dimension) const {
    VectorDiType index = this->relativeIndex(dimension, offset,
                                             dimension2, offset2);

    return _pressurePosition(index, position_dimension);
  }

  ScalarType
  pressurePosition(VectorDiType const& relative_index,
                   int const&          position_dimension) const {
    VectorDiType index = this->relativeIndex(relative_index);

    return _pressurePosition(index, position_dimension);
  }

  VectorDsType&
  velocity() const {
    return _memory->velocity(this->indexValue(Dimensions));
  }

  VectorDsType&
  velocity(int const& dimension, int const& offset) const {
    return _memory->velocity(this->relativeGlobalIndex(dimension, offset));
  }

  VectorDsType&
  velocity(int const& dimension,
           int const& offset,
           int const& dimension2,
           int const& offset2) const {
    return _memory->velocity(
      this->relativeGlobalIndex(dimension, offset,
                                dimension2, offset2));
  }

  VectorDsType&
  velocity(VectorDiType const& index) const {
    return _memory->velocity(this->relativeGlobalIndex(index));
  }

  ScalarType&
  velocity(int const& dimension) const {
    return _memory->velocity(this->globalIndex(), dimension);
  }

  ScalarType&
  velocity(int const& dimension,
           int const& offset,
           int const& velocity_dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(dimension, offset),
                             velocity_dimension);
  }

  ScalarType&
  velocity(int const& dimension,
           int const& offset,
           int const& dimension2,
           int const& offset2,
           int const& velocity_dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(dimension,
                                                       offset,
                                                       dimension2,
                                                       offset2),
                             velocity_dimension);
  }

  ScalarType&
  velocity(VectorDiType const& index,
           int const&          velocity_dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(index),
                             velocity_dimension);
  }

  ScalarType&
  pressure() const {
    return _memory->pressure(this->indexValue(Dimensions));
  }

  ScalarType&
  pressure(int const& dimension, int const& offset) const {
    return _memory->pressure(this->relativeGlobalIndex(dimension, offset));
  }

  ScalarType&
  pressure(int const& dimension,
           int const& offset,
           int const& dimension2,
           int const& offset2) const {
    return _memory->pressure(
      this->relativeGlobalIndex(dimension, offset,
                                dimension2, offset2));
  }

  ScalarType&
  pressure(VectorDiType const& index) const {
    return _memory->pressure(this->relativeGlobalIndex(index));
  }

  int&
  positionInRespectToGeometry() const {
    return _memory->position(this->indexValue(Dimensions));
  }

  int&
  relativePositionInRespectToGeometry(int const& dimension,
                                      int const& offset) const {
    return _memory->position(
      this->relativeGlobalIndex(dimension, offset));
  }

  int&
  positionInRespectToGeometry(int const& dimension,
                              int const& offset,
                              int const& dimension2,
                              int const& offset2) const {
    return _memory->position(
      this->relativeGlobalIndex(dimension, offset,
                                dimension2, offset2));
  }

  int&
  positionInRespectToGeometry(VectorDiType const& index) const {
    return _memory->position(this->relativeGlobalIndex(index));
  }

  VectorDsType&
  fgh() const {
    return _memory->fgh(this->indexValue(Dimensions));
  }

  VectorDsType&
  fgh(int const& dimension, int const& offset) const {
    return _memory->fgh(this->relativeGlobalIndex(dimension, offset));
  }

  VectorDsType&
  fgh(int const& dimension,
      int const& offset,
      int const& dimension2,
      int const& offset2) const {
    return _memory->fgh(
      this->relativeGlobalIndex(dimension, offset,
                                dimension2, offset2));
  }

  VectorDsType&
  fgh(VectorDiType const& index) const {
    return _memory->fgh(this->relativeGlobalIndex(index));
  }

  ScalarType&
  fgh(int const& dimension) const {
    return _memory->fgh(this->globalIndex(), dimension);
  }

  ScalarType&
  fgh(int const& dimension,
      int const& offset,
      int const& fgh_dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(dimension, offset),
                        fgh_dimension);
  }

  ScalarType&
  fgh(int const& dimension,
      int const& offset,
      int const& dimension2,
      int const& offset2,
      int const& fgh_dimension) const {
    return _memory->fgh(
      this->relativeGlobalIndex(dimension, offset,
                                dimension2, offset2),
      fgh_dimension);
  }

  ScalarType&
  fgh(VectorDiType const& index,
      int const&          fgh_dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(index), fgh_dimension);
  }

  VectorDsType&
  convection() const {
    return _memory->convection(this->indexValue(Dimensions));
  }

  ScalarType&
  convection(int const& dimension) const {
    return _memory->convection(this->indexValue(Dimensions), dimension);
  }

  VectorDsType
  absoluteWidth(VectorDiType const& index) const {
    return _memory->gridGeometry()->cellWidth(index);
  }

  ScalarType
  absoluteWidth(VectorDiType const& index,
                int const&          width_dimension) const {
    return _memory->gridGeometry()->cellWidth(index)(width_dimension);
  }

  VectorDsType
  absolutePosition(VectorDiType const& index) const {
    return _memory->gridGeometry()->cellPosition(index);
  }

  ScalarType
  absolutePosition(VectorDiType const& index,
                   int const&          position_dimension) const {
    return _memory->gridGeometry()->cellPosition(index)(position_dimension);
  }

  VectorDsType
  absoluteVelocityPosition(VectorDiType const& absolute_index,
                           int const&          velocity_dimension) const {
    return _velocityPosition(absolute_index,
                             velocity_dimension);
  }

  VectorDsType
  absolutePressurePosition(VectorDiType const& absolute_index) const {
    return _pressurePosition(absolute_index);
  }

  VectorDsType&
  absoluteVelocity(VectorDiType const& index) const {
    return _memory->velocity(this->absoluteGlobalIndex(index));
  }

  ScalarType&
  absoluteVelocity(VectorDiType const& index,
                   int const&          velocity_dimension) const {
    return _memory->velocity(this->absoluteGlobalIndex(index),
                             velocity_dimension);
  }

  ScalarType&
  absolutePressure(VectorDiType const& index) const {
    return _memory->pressure(this->absoluteGlobalIndex(index));
  }

  int&
  absolutePositionInRespectToGeometry(VectorDiType const& index) const {
    return _memory->position(this->absoluteGlobalIndex(index));
  }

  VectorDsType&
  absoluteFgh(VectorDiType const& index) const {
    return _memory->fgh(this->absoluteGlobalIndex(index));
  }

  ScalarType&
  absoluteFgh(VectorDiType const& index,
              int const&          fgh_dimension) const {
    return _memory->fgh(this->absoluteGlobalIndex(index), fgh_dimension);
  }

protected:
  VectorDsType
  _velocityPosition(VectorDiType const& index,
                    int const&          velocity_dimension) const {
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(
      index - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == velocity_dimension) {
        result(d) +=  width(d);
      } else {
        result(d) += 0.5 * width(d);
      }
    }

    return result;
  }

  ScalarType
  _velocityPosition(VectorDiType const& index,
                    int const&          velocity_dimension,
                    int const&          position_dimension) const {
    ScalarType width = _memory->gridGeometry()->cellWidth(index)
                         (position_dimension);
    ScalarType result
      = _memory->gridGeometry()->cellPosition(
      index - _memory->grid()->leftIndent())
          (position_dimension);

    if (position_dimension == velocity_dimension) {
      result += width;
    } else {
      result += 0.5 * width;
    }

    return result;
  }

  VectorDsType
  _pressurePosition(VectorDiType const& index) const {
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(
      index - _memory->grid()->leftIndent())
        + 0.5 * _memory->gridGeometry()->cellWidth(index);

    return result;
  }

  ScalarType
  _pressurePosition(VectorDiType const& index,
                    int const&          position_dimension) const {
    ScalarType result
      = _memory->gridGeometry()->cellPosition(
      index - _memory->grid()->leftIndent())
          (position_dimension)
        + 0.5 * _memory->gridGeometry()->cellWidth(index)(position_dimension);

    return result;
  }

  ScalarType&
  _attribute(int const& global_index,
             int const& attribute_index,
             int const& dimension) const {
    return _memory->attribute(global_index, attribute_index, dimension);
  }

protected:
  MemoryType* _memory;
};
}
}
