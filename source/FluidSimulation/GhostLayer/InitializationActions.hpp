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
  typedef typename TGrid::CellAccessor::VectorDsType Velocity;

  MovingWallVelocityAction(Configuration* parameters,
                           Velocity*      maxVelocity)
    : _configuration(parameters),
    _maxVelocity(maxVelocity) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    for (int d = 0; d < TD; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          current.velocity(d) =
            _configuration->walls[TDimension][TDirection]->velocity() (d);
        } else {
          neighbor.velocity(d) =
            _configuration->walls[TDimension][TDirection]->velocity() (d);
        }
      } else {
        current.velocity(d) =
          2.0 * _configuration->walls[TDimension][TDirection]->velocity() (d)
          - neighbor.velocity(d);
      }
    }
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (current,  *_maxVelocity);
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (neighbor, *_maxVelocity);
  }

  Configuration* _configuration;
  Velocity*      _maxVelocity;
};

template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class InputVelocityAction {
public:
  typedef typename TGrid::CellAccessor::VectorDsType Velocity;

  InputVelocityAction(Configuration* parameters,
                      Velocity*      maxVelocity)
    : _configuration(parameters),
    _maxVelocity(maxVelocity) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    for (int d = 0; d < TD; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          current.velocity(d) =
            _configuration->walls[TDimension][TDirection]->velocity() (d);
        } else {
          neighbor.velocity(d) =
            _configuration->walls[TDimension][TDirection]->velocity() (d);
        }
      } else {
        current.velocity(d) =
          2.0 * _configuration->walls[TDimension][TDirection]->velocity() (d)
          - neighbor.velocity(d);
      }
    }
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (current,  *_maxVelocity);
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (neighbor, *_maxVelocity);
  }

  Configuration* _configuration;
  Velocity*      _maxVelocity;
};

template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class ParabolicInputVelocityAction {
public:
  typedef typename TGrid::CellAccessor::VectorDsType Velocity;

  ParabolicInputVelocityAction(Configuration* parameters,
                               Velocity*      maxVelocity)
    : _configuration(parameters),
    _maxVelocity(maxVelocity) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    for (int d = 0; d < TD; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          current.velocity(d) =
            computeParabolicInputVelocity(
              current,
              _configuration->walls[TDimension][TDirection]->velocity(),
              TDimension);
        } else {
          neighbor.velocity(d) =
            computeParabolicInputVelocity(
              neighbor,
              _configuration->walls[TDimension][TDirection]->velocity(),
              TDimension);
        }
      } else {
        current.velocity(d) =
          2.0 * _configuration->walls[TDimension][TDirection]->velocity() (d)
          - neighbor.velocity(d);
      }
    }
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (current,  *_maxVelocity);
    computeMaxVelocity<typename TGrid::CellAccessor, TScalar, TD>
      (neighbor, *_maxVelocity);
  }

  Configuration* _configuration;
  Velocity*      _maxVelocity;
};

template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class OutputVelocityAction {
public:
 typedef typename TGrid::CellAccessor::VectorDsType Velocity;

  OutputVelocityAction(Velocity* maxVelocity)
    : _maxVelocity(maxVelocity) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    for (int d = 0; d < TD; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          current.velocity(d) =
            neighbor.velocity(d);
        } else {
          neighbor.velocity(d) =
            neighbor.relativeVelocity(d, -1, d);
        }
      } else {
        current.velocity(d) =
          neighbor.velocity(d);
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
  MovingWallFghAction(Configuration* parameters) : _configuration(parameters) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    if (TDirection == 0) {
      current.fgh(TDimension) =
        _configuration->walls[TDimension][TDirection]->velocity() (TDimension);
    } else {
      neighbor.fgh(TDimension) =
        _configuration->walls[TDimension][TDirection]->velocity() (TDimension);
    }
  }

  Configuration* _configuration;
};

template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class InputFghAction {
public:
  InputFghAction(Configuration* parameters) : _configuration(parameters) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    if (TDirection == 0) {
      current.fgh(TDimension) =
        _configuration->walls[TDimension][TDirection]->velocity() (TDimension);
    } else {
      neighbor.fgh(TDimension) =
        _configuration->walls[TDimension][TDirection]->velocity() (TDimension);
    }
  }

  Configuration* _configuration;
};

template <typename TGrid,
          typename TScalar,
          int TD,
          int TDimension,
          int TDirection>
class ParabolicInputFghAction {
public:
  ParabolicInputFghAction(Configuration* parameters)
    : _configuration(parameters) {}

  void
  setValue(typename TGrid::CellAccessor const& current,
           typename TGrid::CellAccessor const& neighbor) {
    if (TDirection == 0) {
      current.fgh(TDimension) =
        computeParabolicInputVelocity(
          current,
          _configuration->walls[TDimension][TDirection]->velocity(),
          TDimension);
    } else {
      neighbor.fgh(TDimension) =
        computeParabolicInputVelocity(
          neighbor,
          _configuration->walls[TDimension][TDirection]->velocity(),
          TDimension);
    }
  }

  Configuration* _configuration;
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
      current.fgh(TDimension) =
        neighbor.fgh(TDimension);
    } else {
      // neighbor.fgh(d) =
      // neighbor.leftCellInDimension(d)->fgh(d);
    }
  }
};
}
}
}
}

#endif
