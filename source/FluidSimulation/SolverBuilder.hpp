#pragma once

#include "Configuration.hpp"
#include "GhostLayer/CornerVelocityHandler.hpp"

#include <functional>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverBuilderTraits>
class SolverBuilder {
private:
  using SolverBuilderTraitsType = TSolverBuilderTraits;

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = typename SolverBuilderTraitsType::template GhostHandlersBuilderType
            <TDimension, TDirection>;

  using SolverTraitsType = typename SolverBuilderTraitsType::SolverTraitsType;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

public:
  using SolverType
          = typename SolverTraitsType::SolverType;

public:
  SolverBuilder(Configuration* configuration,
                SolverType*    solver)
    : _configuration(configuration) {
    _solver = solver;

    VectorDiType processor_size;
    VectorDiType global_cell_size;
    VectorDsType geometry_width;
    processor_size(0)
      = static_cast<ScalarType>(configuration->parallelizationSize(0));
    global_cell_size(0)
      = static_cast<ScalarType>(configuration->size(0));
    geometry_width(0)
      = static_cast<ScalarType>(configuration->width(0));
    processor_size(1)
      = static_cast<ScalarType>(configuration->parallelizationSize(1));
    global_cell_size(1)
      = static_cast<ScalarType>(configuration->size(1));
    geometry_width(1)
      = static_cast<ScalarType>(configuration->width(1));

    if (Dimensions == 3) {
      processor_size(2)
        = static_cast<ScalarType>(configuration->parallelizationSize(2));
      global_cell_size(2)
        = static_cast<ScalarType>(configuration->size(2));
      geometry_width(2)
        = static_cast<ScalarType>(configuration->width(2));
    }

    _solver->memory()->parameters()->re() = configuration->re;

    _solver->memory()->parameters()->gamma() = configuration->gamma;

    _solver->memory()->parameters()->tau() = configuration->tau;

    _solver->memory()->parameters()->alpha() = configuration->alpha;

    _solver->memory()->parameters()->outerLayerSize()
      = configuration->outerLayerSize;

    _solver->memory()->parameters()->innerLayerSize()
      = configuration->innerLayerSize;

    _solver->memory()->parameters()->g(0) = configuration->environment(0);

    _solver->memory()->parameters()->g(1) = configuration->environment(1);

    if (Dimensions == 3) {
      _solver->memory()->parameters()->g(2)
        = configuration->environment(2);
    }

    _solver->memory()->initialize(processor_size,
                                  global_cell_size,
                                  geometry_width);

    setLeftWallAs(configuration->walls[0][0]->type());
    setRightWallAs(configuration->walls[0][1]->type());
    setBottomWallAs(configuration->walls[1][0]->type());
    setTopWallAs(configuration->walls[1][1]->type());
    setBackWallAs(configuration->walls[2][0]->type());
    setFrontWallAs(configuration->walls[2][1]->type());

    auto handlers
      = std::make_shared
        < GhostLayer::CornerVelocityHandlers < SolverTraitsType >> (
      _solver->memory(),
      configuration);

    _solver->ghostHandlers()->cornersHandler
      = std::bind(
      &GhostLayer::CornerVelocityHandlers<SolverTraitsType>::execute,
      handlers);
  }

  void
  setLeftWallAs(WallEnum type) {
    GhostHandlersBuilderType<0, 0> handlers(_configuration, _solver);

    if (_solver->memory()->parallelDistribution()->neighbors[0][0] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setRightWallAs(WallEnum type) {
    GhostHandlersBuilderType<0, 1> handlers(_configuration, _solver);

    if (_solver->memory()->parallelDistribution()->neighbors[0][1] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setBottomWallAs(WallEnum type) {
    if (Dimensions < 2) { return; }

    GhostHandlersBuilderType<1, 0> handlers(_configuration, _solver);

    if (_solver->memory()->parallelDistribution()->neighbors[1][0] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setTopWallAs(WallEnum type) {
    if (Dimensions < 2) { return; }

    GhostHandlersBuilderType<1, 1> handlers(_configuration, _solver);

    if (_solver->memory()->parallelDistribution()->neighbors[1][1] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setBackWallAs(WallEnum type) {
    if (Dimensions < 3) { return; }

    GhostHandlersBuilderType<2, 0> handlers(_configuration, _solver);

    if (_solver->memory()->parallelDistribution()->neighbors[2][0] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setFrontWallAs(WallEnum type) {
    if (Dimensions < 3) { return; }

    GhostHandlersBuilderType<2, 1> handlers(_configuration, _solver);

    if (_solver->memory()->parallelDistribution()->neighbors[2][1] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

private:
  Configuration* _configuration;
  SolverType*    _solver;
};
}
}
