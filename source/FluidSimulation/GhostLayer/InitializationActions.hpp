#ifndef FsiSimulation_FluidSimulation_GhostLayer_Initialization_Actions_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_Initialization_Actions_hpp

#include "FluidSimulation/Configuration.hpp"
#include "FluidSimulation/ParallelDistribution.hpp"
#include "FluidSimulation/functions.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace Initialization {
template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class MovingWallVelocityAction {
public:
  typedef typename TGrid::CellAccessor::CellType::Velocity Velocity;

  MovingWallVelocityAction(Configuration& parameters,
                           Velocity*      maxVelocity)
    : _parameters(parameters),
      _maxVelocity(maxVelocity) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    for (int d = 0; d < TD; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          current.currentCell()->velocity(d) =
            _parameters.walls.velocities[TDimension][TDirection](d);
        } else {
          neighbor.currentCell()->velocity(d) =
            _parameters.walls.velocities[TDimension][TDirection](d);
        }
      } else {
        current.currentCell()->velocity(d) =
          2.0 * _parameters.walls.velocities[TDimension][TDirection](d)
          - neighbor.currentCell()->velocity(d);
      }
    }
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (current,  *_maxVelocity);
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (neighbor, *_maxVelocity);
  }

  Configuration& _parameters;
  Velocity*      _maxVelocity;
};

template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class MovingWallFghAction {
public:
  MovingWallFghAction(Configuration& parameters) : _parameters(parameters) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    if (TDirection == 0) {
      current.currentCell()->fgh(TDimension) =
        _parameters.walls.velocities[TDimension][TDirection](TDimension);
    } else {
      neighbor.currentCell()->fgh(TDimension) =
        _parameters.walls.velocities[TDimension][TDirection](TDimension);
    }
  }

  Configuration& _parameters;
};
}
}
}
}

#endif
