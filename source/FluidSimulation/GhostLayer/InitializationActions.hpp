#ifndef FsiSimulation_FluidSimulation_GhostLayer_Initialization_Actions_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_Initialization_Actions_hpp

#include "Private/utilities.hpp"

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
class InputVelocityAction {
public:
  typedef typename TGrid::CellAccessor::CellType::Velocity Velocity;

  InputVelocityAction(Configuration& parameters,
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
class ParabolicInputVelocityAction {
public:
  typedef typename TGrid::CellAccessor::CellType::Velocity Velocity;

  ParabolicInputVelocityAction(Configuration& parameters,
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
            computeParabolicInputVelocity(
              current,
              _parameters.walls.velocities[TDimension][TDirection],
              TDimension);
        } else {
          neighbor.currentCell()->velocity(d) =
            computeParabolicInputVelocity(
              neighbor,
              _parameters.walls.velocities[TDimension][TDirection],
              TDimension);
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
class OutputVelocityAction {
public:
  typedef typename TGrid::CellAccessor::CellType::Velocity Velocity;

  OutputVelocityAction(Velocity* maxVelocity)
    : _maxVelocity(maxVelocity) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    for (int d = 0; d < TD; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          current.currentCell()->velocity(d) =
            neighbor.currentCell()->velocity(d);
        } else {
          neighbor.currentCell()->velocity(d) =
            neighbor.leftCellInDimension(d)->velocity(d);
        }
      } else {
        current.currentCell()->velocity(d) =
          neighbor.currentCell()->velocity(d);
      }
    }
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (current,  *_maxVelocity);
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (neighbor, *_maxVelocity);
  }

  Velocity* _maxVelocity;
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

template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class InputFghAction {
public:
  InputFghAction(Configuration& parameters) : _parameters(parameters) {}

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

template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class ParabolicInputFghAction {
public:
  ParabolicInputFghAction(Configuration& parameters) : _parameters(
                                                         parameters) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    if (TDirection == 0) {
      current.currentCell()->fgh(TDimension) =
        computeParabolicInputVelocity(
          current,
          _parameters.walls.velocities[TDimension][TDirection],
          TDimension);
    } else {
      neighbor.currentCell()->fgh(TDimension) =
        computeParabolicInputVelocity(
          neighbor,
          _parameters.walls.velocities[TDimension][TDirection],
          TDimension);
    }
  }

  Configuration& _parameters;
};

template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class OutputFghAction {
public:
  OutputFghAction() {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    if (TDirection == 0) {
      current.currentCell()->fgh(TDimension) =
        neighbor.currentCell()->fgh(TDimension);
    } else {
      // neighbor.currentCell()->fgh(d) =
      // neighbor.leftCellInDimension(d)->fgh(d);
    }
  }
};
}
}
}
}

#endif
