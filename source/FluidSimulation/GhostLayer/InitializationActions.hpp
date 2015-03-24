#pragma once

#include "Private/utilities.hpp"

#include "FluidSimulation/Configuration.hpp"
#include "FluidSimulation/ParallelDistribution.hpp"
#include "FluidSimulation/functions.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace Initialization {
template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class MovingWallVelocityAction {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridType = typename SolverTraitsType::GridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  MovingWallVelocityAction(Configuration* parameters,
                           VectorDsType*  maxVelocity)
    : _configuration(parameters),
    _maxVelocity(maxVelocity) {}

  void
  setValue(CellAccessorType const& current,
           CellAccessorType const& neighbor) {
    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          current.velocity(d)
            = _configuration->walls[TDimension][TDirection]->velocity() (d);
        } else {
          neighbor.velocity(d)
            = _configuration->walls[TDimension][TDirection]->velocity() (d);
        }
      } else {
        current.velocity(d)
          = 2.0 * _configuration->walls[TDimension][TDirection]->velocity() (d)
            - neighbor.velocity(d);
      }
    }
    computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
      (current,  *_maxVelocity);
    computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
      (neighbor, *_maxVelocity);
  }

  Configuration* _configuration;
  VectorDsType*  _maxVelocity;
};

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class InputVelocityAction {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridType = typename SolverTraitsType::GridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  InputVelocityAction(Configuration* parameters,
                      VectorDsType*  maxVelocity)
    : _configuration(parameters),
    _maxVelocity(maxVelocity) {}

  void
  setValue(CellAccessorType const& current,
           CellAccessorType const& neighbor) {
    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          current.velocity(d)
            = _configuration->walls[TDimension][TDirection]->velocity() (d);
        } else {
          neighbor.velocity(d)
            = _configuration->walls[TDimension][TDirection]->velocity() (d);
        }
      } else {
        current.velocity(d)
          = 2.0 * _configuration->walls[TDimension][TDirection]->velocity() (d)
            - neighbor.velocity(d);
      }
    }
    computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
      (current,  *_maxVelocity);
    computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
      (neighbor, *_maxVelocity);
  }

  Configuration* _configuration;
  VectorDsType   _;
  VectorDsType*  _maxVelocity;
};

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class ParabolicInputVelocityAction {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridType = typename SolverTraitsType::GridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  ParabolicInputVelocityAction(Configuration* parameters,
                               VectorDsType*  maxVelocity)
    : _configuration(parameters),
    _maxVelocity(maxVelocity) {}

  void
  setValue(CellAccessorType const& current,
           CellAccessorType const& neighbor) {
    ScalarType result;

    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        result
          = _configuration->walls[TDimension][TDirection]->velocity()(
          TDimension);

        if (TDirection == 0) {
          compute_parabolic_input_velocity(current, TDimension, result);
          current.velocity(TDimension) = result;
        } else {
          compute_parabolic_input_velocity(neighbor, TDimension, result);
          neighbor.velocity(TDimension) = result;
        }
      } else {
        current.velocity(d)
          = 2.0 * _configuration->walls[TDimension][TDirection]->velocity() (d)
            - neighbor.velocity(d);
      }
    }
    computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
      (current,  *_maxVelocity);
    computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
      (neighbor, *_maxVelocity);
  }

  Configuration* _configuration;
  VectorDsType*  _maxVelocity;
};

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class OutputVelocityAction {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridType = typename SolverTraitsType::GridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  OutputVelocityAction(VectorDsType* maxVelocity)
    : _maxVelocity(maxVelocity) {}

  void
  setValue(CellAccessorType const& current,
           CellAccessorType const& neighbor) {
    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          current.velocity(d)
            = neighbor.velocity(d);
        } else {
          neighbor.velocity(d)
            = neighbor.velocity(d, -1, d);
        }
      } else {
        current.velocity(d)
          = neighbor.velocity(d);
      }
    }
    computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
      (current,  *_maxVelocity);
    computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
      (neighbor, *_maxVelocity);
  }

  VectorDsType* _maxVelocity;
};

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class MovingWallFghAction {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridType = typename SolverTraitsType::GridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

public:
  MovingWallFghAction(Configuration* parameters) : _configuration(parameters) {}

  void
  setValue(CellAccessorType const& current,
           CellAccessorType const& neighbor) {
    if (TDirection == 0) {
      current.fgh(TDimension)
        = _configuration->walls[TDimension][TDirection]->velocity() (
        TDimension);
    } else {
      neighbor.fgh(TDimension)
        = _configuration->walls[TDimension][TDirection]->velocity() (
        TDimension);
    }
  }

  Configuration* _configuration;
};

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class InputFghAction {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridType = typename SolverTraitsType::GridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  InputFghAction(Configuration* parameters) : _configuration(parameters) {}

  void
  setValue(CellAccessorType const& current,
           CellAccessorType const& neighbor) {
    if (TDirection == 0) {
      current.fgh(TDimension)
        = _configuration->walls[TDimension][TDirection]->velocity() (
        TDimension);
    } else {
      neighbor.fgh(TDimension)
        = _configuration->walls[TDimension][TDirection]->velocity() (
        TDimension);
    }
  }

  Configuration* _configuration;
};

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class ParabolicInputFghAction {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridType = typename SolverTraitsType::GridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  ParabolicInputFghAction(Configuration* parameters)
    : _configuration(parameters) {}

  void
  setValue(CellAccessorType const& current,
           CellAccessorType const& neighbor) {
    ScalarType result;
    result
      = _configuration->walls[TDimension][TDirection]->velocity()(TDimension);

    if (TDirection == 0) {
      compute_parabolic_input_velocity(current, TDimension, result);
      current.fgh(TDimension) = result;
    } else {
      compute_parabolic_input_velocity(neighbor, TDimension, result);
      neighbor.fgh(TDimension) = result;
    }
  }

  Configuration* _configuration;
};

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class OutputFghAction {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridType = typename SolverTraitsType::GridType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

public:
  OutputFghAction() {}

  void
  setValue(CellAccessorType const& current,
           CellAccessorType const& neighbor) {
    if (TDirection == 0) {
      current.fgh(TDimension) = neighbor.fgh(TDimension);
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
