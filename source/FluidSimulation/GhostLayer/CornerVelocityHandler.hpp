#pragma once

#include "Private/utilities.hpp"

#include "FluidSimulation/Configuration.hpp"

#include <Uni/ExecutionControl/exception>
#include <Uni/Logging/macros>

#include <functional>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
/**
 * direction --- direction (0 or 1) per dimension,
 *               describe the corner location.
 * index --- index of the root cell of that offset is made.
 * offset --- offset to make from the root cell index.
 */
template <typename TSolverTraits,
          int TDimension1,
          int TDirection1,
          int TDimension2,
          int TDirection2>
class CornerVelocityHandler {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using MemoryType = typename SolverTraitsType::MemoryType;

  using GridType = typename SolverTraitsType::GridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

public:
  CornerVelocityHandler(MemoryType const*    memory,
                        Configuration const* configuration) :
    _memory(memory),
    _configuration(configuration),
    _dimension(-1),
    _position(VectorDsType::Zero()) {
    for (int d = 0; d < Dimensions; ++d) {
      if (d != TDimension1
          && d != TDimension2) {
        _dimension = d;
        break;
      }
    }
    composeIndexAndOffset(TDimension1, TDirection1);
    composeIndexAndOffset(TDimension2, TDirection2);
  }

  VectorDsType
  getVelocity(int const&              dimension,
              int const&              direction,
              int const&              dimension2,
              int const&              direction2,
              VectorDsType const&     position,
              CellAccessorType const& accessor) const {
    auto result = VectorDsType::Zero().eval();

    if ((_configuration->walls[dimension][direction]->type()
         == WallEnum::Input)) {
      result(dimension)
        = 2 * _configuration->walls[dimension][direction]
          ->velocity()(dimension);

      result(dimension2)
        = 2 * _configuration->walls[dimension][direction]
          ->velocity()(dimension2);
    } else if ((_configuration->walls[dimension][direction]->type()
                == WallEnum::ParabolicInput)) {
      result(dimension)
        =  _configuration->walls[dimension][direction]->velocity()(dimension);
      compute_parabolic_input(
        position,
        _memory->gridGeometry()->size(),
        dimension,
        result(dimension));
        result(dimension) = 2 * result(dimension);

        result(dimension2)
        = 2 * _configuration->walls[dimension][direction]
          ->velocity()(dimension2);

      return result;
    } else if ((_configuration->walls[dimension][direction]->type()
                == WallEnum::Output)) {
      result(dimension)
        = accessor.velocity(dimension,
                            _offset(dimension) + _offset2(dimension),
                            dimension);
      result(dimension)
        += accessor.velocity(dimension,
                             _offset(dimension) + _offset2(dimension),
                             dimension2,
                             _offset2(dimension2),
                             dimension);
      result(dimension2)
        = accessor.velocity(dimension2,
                            _offset(dimension2),
                            dimension,
                            _offset2(dimension),
                            dimension2);
      result(dimension2)
        = accessor.velocity(dimension2,
                            _offset(dimension2),
                            dimension,
                            2 * _offset2(dimension),
                            dimension2);
    }
    result(dimension)
      -= accessor.velocity(dimension,
                           _offset(dimension),
                           dimension2,
                           _offset2(dimension2),
                           dimension);
    result(dimension2)
      -= accessor.velocity(dimension2,
                           _offset(dimension2),
                           dimension,
                           _offset2(dimension),
                           dimension2);

    return result;
  }

  void
  computeCornerVelocity(CellAccessorType const& accessor,
                        VectorDsType const&     position) const {
    VectorDsType corner_velocity;

    corner_velocity = getVelocity(TDimension1,
                                  TDirection1,
                                  TDimension2,
                                  TDirection2,
                                  position,
                                  accessor);
    corner_velocity += getVelocity(TDimension2,
                                   TDirection2,
                                   TDimension1,
                                   TDirection1,
                                   position,
                                   accessor);
    corner_velocity /= 2.0;
    accessor.velocity(TDimension1, _offset(TDimension1), TDimension1)
      = corner_velocity(TDimension1);
    accessor.velocity(TDimension2, _offset(TDimension2), TDimension2)
      = corner_velocity(TDimension2);
  }

  void
  computeCornerVelocity() const {
    auto it = _memory->grid()->innerGrid.at(_index);

    if (_dimension == -1) {
      computeCornerVelocity(*it, _position);
    } else {
      auto position = _position;

      for (; it->indexValue(_dimension)
           < _memory->grid()->innerGrid.end()->indexValue(_dimension);
           it->index(_dimension) += 1) {
        position(_dimension) = it->position(_dimension);
        computeCornerVelocity(*it, position);
      }
    }
  }

private:
  void
  _computeAverageVelocity(CellAccessorType const& accessor) {
    auto result = VectorDsType::Zero();
    auto offset = VectorDiType::Zero();
    int  d      = 0;

    while (true) {
      if (offset(d) == 0) {
        ++offset(d);
      } else {
        ++d;

        if (d == Dimensions) {
          break;
        }

        for (int i = 0; i < d; ++i) {
          offset(i) = 0;
        }
        continue;
      }

      auto result2 = VectorDsType::Zero();

      for (int d = 0; d < Dimensions; ++d) {
        result2(d) = accessor.velocity(offset, d);
        auto offset2 = offset;
        offset2(d) -= 1;
        result2(d) += accessor.velocity(offset2, d);
        result2(d) *= 0.5;
      }
      result += result2;
    }

    if (Dimensions == 2) {
      result = result / 4.0;
    } else if (Dimensions == 3) {
      result = result / 8.0;
    }

    return result;
  }

  void
  composeIndexAndOffset(int const& dimension,
                        int const& direction) {
    if (direction == 0) {
      _index(dimension)
        = _memory->grid()->innerGrid.leftIndent(dimension) -
          1;
      _offset(dimension)  = 0;
      _offset2(dimension) = +1;
    } else if (direction == 1) {
      _index(dimension)
                           = _memory->grid()->innerLimit(dimension) - 1;
      _offset(dimension)   = -1;
      _offset2(dimension)  = -1;
      _position(dimension) = _memory->gridGeometry()->size(dimension);
    } else {
      throwException("Unknown direction value");
    }
  }
  MemoryType const*    _memory;
  Configuration const* _configuration;
  int                  _dimension;

  VectorDiType _index;
  VectorDiType _offset;
  VectorDiType _offset2;
  VectorDsType _position;
};

template <typename TSolverTraits,
          int TDimensions = TSolverTraits::Dimensions>
struct CornerVelocityHandlers {};

template <typename TSolverTraits>
struct CornerVelocityHandlers<TSolverTraits, 2> {
  using SolverTraitsType = TSolverTraits;

  using MemoryType = typename SolverTraitsType::MemoryType;

  CornerVelocityHandlers(MemoryType const*    memory,
                         Configuration const* configuration) :
    _rb(memory, configuration),
    _lt(memory, configuration) {}

  void
  execute() const {
    _rb.computeCornerVelocity();
    _lt.computeCornerVelocity();
  }

private:
  CornerVelocityHandler<TSolverTraits, 0, 1, 1, 0> _rb;
  CornerVelocityHandler<TSolverTraits, 0, 0, 1, 1> _lt;
};

template <typename TSolverTraits>
struct CornerVelocityHandlers<TSolverTraits, 3> {};
}
}
}
