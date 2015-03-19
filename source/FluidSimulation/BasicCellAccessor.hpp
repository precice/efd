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
  attribute(short const&    attribute_index,
            unsigned const& dimension = 0) const {
    return _attribute(this->globalIndex(), attribute_index, dimension);
  }

  VectorDsType
  width() const {
    return _memory->gridGeometry()->cellWidth(this->indexValues()
                                              - _memory->grid()->leftIndent());
  }

  ScalarType
  width(int const& dimension) const {
    return _memory->gridGeometry()->cellWidth(
      this->indexValues() - _memory->grid()->leftIndent())(dimension);
  }

  VectorDsType
  position() const {
    return _memory->gridGeometry()->cellPosition(this->indexValues()
                                                 - _memory->grid()->leftIndent());
  }

  ScalarType
  position(int const& dimension) const {
    return _memory->gridGeometry()->cellPosition(
      this->indexValues()
      - _memory->grid()->leftIndent())(dimension);
  }

  VectorDsType
  velocityPosition(int const& dimension) const {
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(this->indexValues()
                                              - _memory->grid()->leftIndent())
        + 0.5 * _memory->gridGeometry()->cellWidth(this->indexValues());
    result(dimension) +=
      0.5 * _memory->gridGeometry()->cellWidth(
        this->indexValues()) (dimension);

    return result;
  }

  ScalarType
  velocityPosition(int const& dimension,
                   int const& position_dimension) const {
    ScalarType result
      = _memory->gridGeometry()->cellPosition(
      this->indexValues() - _memory->grid()->leftIndent())
          (position_dimension);

    if (position_dimension == dimension) {
      result +=
        _memory->gridGeometry()->cellWidth(this->indexValues())
          (position_dimension);
    } else {
      result +=
        0.5 * _memory->gridGeometry()->cellWidth(this->indexValues())
          (position_dimension);
    }

    return result;
  }

  VectorDsType
  pressurePosition() const {
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(this->indexValues()
                                              - _memory->grid()->leftIndent())
        + 0.5 * _memory->gridGeometry()->cellWidth(this->indexValues());

    return result;
  }

  VectorDsType&
  velocity() const {
    return _memory->velocity(this->indexValue(Dimensions));
  }

  ScalarType&
  velocity(int const& dimension) const {
    return _memory->velocity(this->indexValue(Dimensions), dimension);
  }

  ScalarType&
  pressure() const {
    return _memory->pressure(this->indexValue(Dimensions));
  }

  int&
  positionInRespectToGeometry() const {
    return _memory->position(this->indexValue(Dimensions));
  }

  VectorDsType&
  fgh() const {
    return _memory->fgh(this->indexValue(Dimensions));
  }

  ScalarType&
  fgh(int const& dimension) const {
    return _memory->fgh(this->indexValue(Dimensions), dimension);
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
    VectorDsType width = _memory->gridGeometry()->cellWidth(absolute_index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(absolute_index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == velocity_dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  absolutePressurePosition(VectorDiType const& absolute_index) const {
    VectorDsType width = _memory->gridGeometry()->cellWidth(absolute_index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(absolute_index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
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

  VectorDsType
  relativeWidth(int const& dimension, int const& offset) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(dimension, offset)
      - _memory->grid()->leftIndent());
  }

  ScalarType
  relativeWidth(int const& dimension,
                int const& offset,
                int const& width_dimension) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(dimension, offset)
      - _memory->grid()->leftIndent())(width_dimension);
  }

  VectorDsType
  relativePosition(int const& dimension, int const& offset) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(dimension, offset)
      - _memory->grid()->leftIndent());
  }

  ScalarType
  relativePosition(int const& dimension,
                   int const& offset,
                   int const& position_dimension) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(dimension, offset)
      - _memory->grid()->leftIndent())(position_dimension);
  }

  VectorDsType
  relativeVelocityPosition(int const& dimension,
                           int const& offset,
                           int const& velocity_dimension) const {
    VectorDiType index = this->relativeIndex(dimension, offset);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == velocity_dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  relativePressurePosition(int const& dimension, int const& offset) const {
    VectorDiType index = this->relativeIndex(dimension, offset);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  relativeVelocity(int const& dimension, int const& offset) const {
    return _memory->velocity(this->relativeGlobalIndex(dimension, offset));
  }

  ScalarType&
  relativeVelocity(int const& dimension,
                   int const& offset,
                   int const& velocity_dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(dimension, offset),
                             velocity_dimension);
  }

  ScalarType&
  relativePressure(int const& dimension, int const& offset) const {
    return _memory->pressure(this->relativeGlobalIndex(dimension, offset));
  }

  int&
  relativePositionInRespectToGeometry(int const& dimension, int
                                      const& offset) const {
    return _memory->position(this->relativeGlobalIndex(
                               dimension, offset));
  }

  VectorDsType&
  relativeFgh(int const& dimension, int const& offset) const {
    return _memory->fgh(this->relativeGlobalIndex(dimension, offset));
  }

  ScalarType&
  relativeFgh(int const& dimension,
              int const& offset,
              int const& fgh_dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(dimension, offset),
                        fgh_dimension);
  }

  VectorDsType
  relativeWidth(int const& dimension,
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

  ScalarType
  relativeWidth(int const& dimension,
                int const& offset,
                int const& dimension2,
                int const& offset2,
                int const& width_dimension) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(dimension,
                          offset,
                          dimension2,
                          offset2)
      - _memory->grid()->leftIndent())(width_dimension);
  }

  VectorDsType
  relativePosition(int const& dimension,
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

  ScalarType
  relativePosition(int const& dimension,
                   int const& offset,
                   int const& dimension2,
                   int const& offset2,
                   int const& position_dimension) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(dimension,
                          offset,
                          dimension2,
                          offset2)
      - _memory->grid()->leftIndent())(position_dimension);
  }

  VectorDsType
  relativeVelocityPosition(int const& dimension,
                           int const& offset,
                           int const& dimension2,
                           int const& offset2,
                           int const& velocity_dimension) const {
    VectorDiType index = this->relativeIndex(dimension,
                                             offset,
                                             dimension2,
                                             offset2);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == velocity_dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  relativePressurePosition(int const& dimension,
                           int const& offset,
                           int const& dimension2,
                           int const& offset2) const {
    VectorDiType index = this->relativeIndex(dimension,
                                             offset,
                                             dimension2,
                                             offset2);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  relativeVelocity(int const& dimension,
                   int const& offset,
                   int const& dimension2,
                   int const& offset2) const {
    return _memory->velocity(this->relativeGlobalIndex(dimension,
                                                       offset,
                                                       dimension2,
                                                       offset2));
  }

  ScalarType&
  relativeVelocity(int const& dimension,
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
  relativePressure(int const& dimension,
                   int const& offset,
                   int const& dimension2,
                   int const& offset2) const {
    return _memory->pressure(this->relativeGlobalIndex(dimension,
                                                       offset,
                                                       dimension2,
                                                       offset2));
  }

  int&
  relativePositionInRespectToGeometry(int const& dimension,
                                      int const& offset,
                                      int const& dimension2,
                                      int const& offset2) const {
    return _memory->position(this->relativeGlobalIndex(
                               dimension,
                               offset,
                               dimension2,
                               offset2));
  }

  VectorDsType&
  relativeFgh(int const& dimension, int const& offset, int const& dimension2,
              int const& offset2) const {
    return _memory->fgh(this->relativeGlobalIndex(dimension,
                                                  offset,
                                                  dimension2,
                                                  offset2));
  }

  ScalarType&
  relativeFgh(int const& dimension, int const& offset, int const& dimension2,
              int const& offset2,
              int const& fgh_dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(dimension,
                                                  offset,
                                                  dimension2,
                                                  offset2),
                        fgh_dimension);
  }

  VectorDsType
  relativeWidth(VectorDiType const& index) const {
    return _memory->gridGeometry()->cellWidth(
      this->relativeIndex(index)
      - _memory->grid()->leftIndent());
  }

  VectorDsType
  relativePosition(VectorDiType const& index) const {
    return _memory->gridGeometry()->cellPosition(
      this->relativeIndex(index)
      - _memory->grid()->leftIndent());
  }

  VectorDsType
  relativeVelocityPosition(VectorDiType const& relative_index,
                           int const&          velocity_dimension) const {
    VectorDiType index = this->relativeIndex(relative_index);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == velocity_dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  relativePressurePosition(VectorDiType const& relative_index) const {
    VectorDiType index = this->relativeIndex(relative_index);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  relativeVelocity(VectorDiType const& index) const {
    return _memory->velocity(this->relativeGlobalIndex(index));
  }

  ScalarType&
  relativeVelocity(VectorDiType const& index,
                   int const&          velocity_dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(index),
                             velocity_dimension);
  }

  ScalarType&
  relativePressure(VectorDiType const& index) const {
    return _memory->pressure(this->relativeGlobalIndex(index));
  }

  int&
  relativePositionInRespectToGeometry(VectorDiType const& index) const {
    return _memory->position(this->relativeGlobalIndex(
                               index));
  }

  VectorDsType&
  relativeFgh(VectorDiType const& index) const {
    return _memory->fgh(this->relativeGlobalIndex(index));
  }

  ScalarType&
  relativeFgh(VectorDiType const& index,
              int const&          fgh_dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(index),
                        fgh_dimension);
  }

  VectorDsType
  leftWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, -1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, -1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  leftPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  leftVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1));
  }

  ScalarType&
  leftVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1), dimension);
  }

  ScalarType&
  leftPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, -1));
  }

  int&
  leftPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0,
                                                       -1));
  }

  VectorDsType&
  leftFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1));
  }

  ScalarType&
  leftFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1), dimension);
  }

  VectorDsType
  rightWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, 1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, 1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  rightPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  rightVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1));
  }

  ScalarType&
  rightVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1), dimension);
  }

  ScalarType&
  rightPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, 1));
  }

  int&
  rightPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0,
                                                       1));
  }

  VectorDsType&
  rightFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1));
  }

  ScalarType&
  rightFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1), dimension);
  }

  VectorDsType
  bottomWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(1, -1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  bottomPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(1, -1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  bottomVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(1, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  bottomPressurePosition() const {
    VectorDiType index = this->relativeIndex(1, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  bottomVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(1, -1));
  }

  ScalarType&
  bottomVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(1, -1), dimension);
  }

  ScalarType&
  bottomPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(1, -1));
  }

  int&
  bottomPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(1,
                                                       -1));
  }

  VectorDsType&
  bottomFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(1, -1));
  }

  ScalarType&
  bottomFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(1, -1), dimension);
  }

  VectorDsType
  topWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(1, 1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  topPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(1, 1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  topVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(1, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  topPressurePosition() const {
    VectorDiType index = this->relativeIndex(1, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  topVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(1, 1));
  }

  ScalarType&
  topVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(1, 1), dimension);
  }

  ScalarType&
  topPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(1, 1));
  }

  int&
  topPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(1,
                                                       1));
  }

  VectorDsType&
  topFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(1, 1));
  }

  ScalarType&
  topFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(1, 1), dimension);
  }

  VectorDsType
  backWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(2, -1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  backPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(2, -1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  backVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  backPressurePosition() const {
    VectorDiType index = this->relativeIndex(2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  backVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(2, -1));
  }

  ScalarType&
  backVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(2, -1), dimension);
  }

  ScalarType&
  backPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(2, -1));
  }

  int&
  backPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(2,
                                                       -1));
  }

  VectorDsType&
  backFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(2, -1));
  }

  ScalarType&
  backFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(2, -1), dimension);
  }

  VectorDsType
  frontWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(2, 1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  frontPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(2, 1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  frontVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  frontPressurePosition() const {
    VectorDiType index = this->relativeIndex(2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  frontVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(2, 1));
  }

  ScalarType&
  frontVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(2, 1), dimension);
  }

  ScalarType&
  frontPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(2, 1));
  }

  int&
  frontPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(2,
                                                       1));
  }

  VectorDsType&
  frontFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(2, 1));
  }

  ScalarType&
  frontFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(2, 1), dimension);
  }

  VectorDsType
  leftBottomWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, -1, 1, -1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftBottomPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, -1, 1,
                                                                     -1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftBottomVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, -1, 1, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  leftBottomPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, -1, 1, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  leftBottomVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1, 1, -1));
  }

  ScalarType&
  leftBottomVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1, 1, -1),
                             dimension);
  }

  ScalarType&
  leftBottomPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, -1, 1, -1));
  }

  int&
  leftBottomPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0, -1,
                                                       1,
                                                       -1));
  }

  VectorDsType&
  leftBottomFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1, 1, -1));
  }

  ScalarType&
  leftBottomFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1, 1, -1), dimension);
  }

  VectorDsType
  leftTopWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, -1, 1, 1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftTopPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, -1, 1,
                                                                     1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftTopVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, -1, 1, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  leftTopPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, -1, 1, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  leftTopVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1, 1, 1));
  }

  ScalarType&
  leftTopVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1, 1, 1), dimension);
  }

  ScalarType&
  leftTopPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, -1, 1, 1));
  }

  int&
  leftTopPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0, -1,
                                                       1,
                                                       1));
  }

  VectorDsType&
  leftTopFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1, 1, 1));
  }

  ScalarType&
  leftTopFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1, 1, 1), dimension);
  }

  VectorDsType
  rightBottomWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, 1, 1, -1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightBottomPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, 1, 1,
                                                                     -1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightBottomVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, 1, 1, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  rightBottomPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, 1, 1, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  rightBottomVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1, 1, -1));
  }

  ScalarType&
  rightBottomVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1, 1, -1), dimension);
  }

  ScalarType&
  rightBottomPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, 1, 1, -1));
  }

  int&
  rightBottomPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0, 1,
                                                       1,
                                                       -1));
  }

  VectorDsType&
  rightBottomFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1, 1, -1));
  }

  ScalarType&
  rightBottomFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1, 1, -1), dimension);
  }

  VectorDsType
  rightTopWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, 1, 1, 1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightTopPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, 1, 1, 1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightTopVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, 1, 1, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  rightTopPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, 1, 1, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  rightTopVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1, 1, 1));
  }

  ScalarType&
  rightTopVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1, 1, 1), dimension);
  }

  ScalarType&
  rightTopPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, 1, 1, 1));
  }

  int&
  rightTopPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0, 1,
                                                       1,
                                                       1));
  }

  VectorDsType&
  rightTopFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1, 1, 1));
  }

  ScalarType&
  rightTopFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1, 1, 1), dimension);
  }

  VectorDsType
  leftBackWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, -1, 2, -1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftBackPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, -1, 2,
                                                                     -1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftBackVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, -1, 2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  leftBackPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, -1, 2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  leftBackVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1, 2, -1));
  }

  ScalarType&
  leftBackVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1, 2, -1),
                             dimension);
  }

  ScalarType&
  leftBackPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, -1, 2, -1));
  }

  int&
  leftBackPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0, -1,
                                                       2,
                                                       -1));
  }

  VectorDsType&
  leftBackFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1, 2, -1));
  }

  ScalarType&
  leftBackFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1, 2, -1), dimension);
  }

  VectorDsType
  leftFrontWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, -1, 2, 1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftFrontPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, -1, 2,
                                                                     1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  leftFrontVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, -1, 2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  leftFrontPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, -1, 2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  leftFrontVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1, 2, 1));
  }

  ScalarType&
  leftFrontVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, -1, 2, 1), dimension);
  }

  ScalarType&
  leftFrontPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, -1, 2, 1));
  }

  int&
  leftFrontPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0, -1,
                                                       2,
                                                       1));
  }

  VectorDsType&
  leftFrontFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1, 2, 1));
  }

  ScalarType&
  leftFrontFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, -1, 2, 1), dimension);
  }

  VectorDsType
  rightBackWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, 1, 2, -1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightBackPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, 1, 2,
                                                                     -1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightBackVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, 1, 2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  rightBackPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, 1, 2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  rightBackVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1, 2, -1));
  }

  ScalarType&
  rightBackVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1, 2, -1), dimension);
  }

  ScalarType&
  rightBackPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, 1, 2, -1));
  }

  int&
  rightBackPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0, 1,
                                                       2,
                                                       -1));
  }

  VectorDsType&
  rightBackFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1, 2, -1));
  }

  ScalarType&
  rightBackFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1, 2, -1), dimension);
  }

  VectorDsType
  rightFrontWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(0, 1, 2, 1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightFrontPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(0, 1, 2, 1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  rightFrontVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(0, 1, 2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  rightFrontPressurePosition() const {
    VectorDiType index = this->relativeIndex(0, 1, 2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  rightFrontVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1, 2, 1));
  }

  ScalarType&
  rightFrontVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(0, 1, 2, 1), dimension);
  }

  ScalarType&
  rightFrontPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(0, 1, 2, 1));
  }

  int&
  rightFrontPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(0, 1,
                                                       2,
                                                       1));
  }

  VectorDsType&
  rightFrontFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1, 2, 1));
  }

  ScalarType&
  rightFrontFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(0, 1, 2, 1), dimension);
  }

  VectorDsType
  bottomBackWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(1, -1, 2, -1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  bottomBackPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(1, -1, 2,
                                                                     -1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  bottomBackVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(1, -1, 2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  bottomBackPressurePosition() const {
    VectorDiType index = this->relativeIndex(1, -1, 2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  bottomBackVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(1, -1, 2, -1));
  }

  ScalarType&
  bottomBackVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(1, -1, 2, -1),
                             dimension);
  }

  ScalarType&
  bottomBackPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(1, -1, 2, -1));
  }

  int&
  bottomBackPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(1, -1,
                                                       2,
                                                       -1));
  }

  VectorDsType&
  bottomBackFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(1, -1, 2, -1));
  }

  ScalarType&
  bottomBackFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(1, -1, 2, -1), dimension);
  }

  VectorDsType
  bottomFrontWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(1, -1, 2, 1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  bottomFrontPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(1, -1, 2,
                                                                     1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  bottomFrontVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(1, -1, 2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  bottomFrontPressurePosition() const {
    VectorDiType index = this->relativeIndex(1, -1, 2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  bottomFrontVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(1, -1, 2, 1));
  }

  ScalarType&
  bottomFrontVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(1, -1, 2, 1), dimension);
  }

  ScalarType&
  bottomFrontPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(1, -1, 2, 1));
  }

  int&
  bottomFrontPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(1, -1,
                                                       2,
                                                       1));
  }

  VectorDsType&
  bottomFrontFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(1, -1, 2, 1));
  }

  ScalarType&
  bottomFrontFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(1, -1, 2, 1), dimension);
  }

  VectorDsType
  topBackWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(1, 1, 2, -1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  topBackPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(1, 1, 2,
                                                                     -1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  topBackVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(1, 1, 2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  topBackPressurePosition() const {
    VectorDiType index = this->relativeIndex(1, 1, 2, -1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  topBackVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(1, 1, 2, -1));
  }

  ScalarType&
  topBackVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(1, 1, 2, -1), dimension);
  }

  ScalarType&
  topBackPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(1, 1, 2, -1));
  }

  int&
  topBackPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(1, 1,
                                                       2,
                                                       -1));
  }

  VectorDsType&
  topBackFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(1, 1, 2, -1));
  }

  ScalarType&
  topBackFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(1, 1, 2, -1), dimension);
  }

  VectorDsType
  topFrontWidth() const {
    return _memory->gridGeometry()->cellWidth(this->relativeIndex(1, 1, 2, 1)
                                              - _memory->grid()->leftIndent());
  }

  VectorDsType
  topFrontPosition() const {
    return _memory->gridGeometry()->cellPosition(this->relativeIndex(1, 1, 2, 1)
                                                 - _memory->grid()->leftIndent());
  }

  VectorDsType
  topFrontVelocityPosition(int const& dimension) const {
    VectorDiType index = this->relativeIndex(1, 1, 2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent());

    for (int d = 0; d < Dimensions; ++d) {
      if (d == dimension) {
        result(d) +=  width;
      } else {
        result(d) += 0.5 * width;
      }
    }

    return result;
  }

  VectorDsType
  topFrontPressurePosition() const {
    VectorDiType index = this->relativeIndex(1, 1, 2, 1);
    VectorDsType width = _memory->gridGeometry()->cellWidth(index);
    VectorDsType result
      = _memory->gridGeometry()->cellPosition(index
                                              - _memory->grid()->leftIndent())
        + 0.5 * width;

    return result;
  }

  VectorDsType&
  topFrontVelocity() const {
    return _memory->velocity(this->relativeGlobalIndex(1, 1, 2, 1));
  }

  ScalarType&
  topFrontVelocity(int const& dimension) const {
    return _memory->velocity(this->relativeGlobalIndex(1, 1, 2, 1), dimension);
  }

  ScalarType&
  topFrontPressure() const {
    return _memory->pressure(this->relativeGlobalIndex(1, 1, 2, 1));
  }

  int&
  topFrontPositionInRespectToGeometry() const {
    return _memory->position(this->relativeGlobalIndex(1, 1, 2, 1));
  }

  VectorDsType&
  topFrontFgh() const {
    return _memory->fgh(this->relativeGlobalIndex(1, 1, 2, 1));
  }

  ScalarType&
  topFrontFgh(int const& dimension) const {
    return _memory->fgh(this->relativeGlobalIndex(1, 1, 2, 1), dimension);
  }

protected:
  ScalarType&
  _attribute(std::size_t const& global_index,
             short const&       attribute_index,
             unsigned const&    dimension) const {
    return _memory->attribute(global_index, attribute_index, dimension);
  }

protected:
  MemoryType* _memory;
};
}
}
