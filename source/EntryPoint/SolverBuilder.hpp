#pragma once

#include "FluidSimulation/SolverBuilderTraits.hpp"

namespace FsiSimulation {
namespace EntryPoint {
template <typename TScalar,
          int TDimensions,
          int TSolverType   = 0,
          int TOutputWriter = 0>
class SolverBuilder {
  static_assert((TDimensions > 1) && (TDimensions < 4),
                "Only 2D and 3D simulations supported");

private:
  using SolverBuilderTraitsType
          = FluidSimulation::SolverBuilderTraits
            <TScalar, TDimensions, TSolverType, TOutputWriter>;

  enum {
    Dimensions = SolverBuilderTraitsType::Dimensions
  };

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = typename SolverBuilderTraitsType::template GhostHandlersBuilderType
            <TDimension, TDirection>;

  using VectorDsType = typename SolverBuilderTraitsType::VectorDsType;

  using VectorDiType = typename SolverBuilderTraitsType::VectorDiType;

public:
  using SolverType = typename SolverBuilderTraitsType::SolverType;

public:
  SolverBuilder(FluidSimulation::Configuration* configuration)
    : _configuration(configuration) {
    _simulation = new SolverType();

    VectorDiType processor_size;
    VectorDiType global_cell_size;
    VectorDsType geometry_width;
    processor_size(0)
      = static_cast<TScalar>(configuration->parallelizationSize(0));
    global_cell_size(0)
      = static_cast<TScalar>(configuration->size(0));
    geometry_width(0)
      = static_cast<TScalar>(configuration->width(0));
    processor_size(1)
      = static_cast<TScalar>(configuration->parallelizationSize(1));
    global_cell_size(1)
      = static_cast<TScalar>(configuration->size(1));
    geometry_width(1)
      = static_cast<TScalar>(configuration->width(1));

    if (Dimensions == 3) {
      processor_size(2)
        = static_cast<TScalar>(configuration->parallelizationSize(2));
      global_cell_size(2)
        = static_cast<TScalar>(configuration->size(2));
      geometry_width(2)
        = static_cast<TScalar>(configuration->width(2));
    }

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    _simulation->_memory.parameters()->re() = configuration->re;

    _simulation->_memory.parameters()->gamma() = configuration->gamma;

    _simulation->_memory.parameters()->tau() = configuration->tau;

    _simulation->_memory.parameters()->alpha() = configuration->alpha;

    _simulation->_memory.parameters()->immersedBoundaryMethod()
      = configuration->immersedBoundaryMethod;

    _simulation->_memory.parameters()->g(0) = configuration->environment(
      0);

    _simulation->_memory.parameters()->g(1) = configuration->environment(
      1);

    _simulation->_iterationLimit = configuration->iterationLimit;
    _simulation->_timeLimit      = configuration->timeLimit;
    _simulation->_plotInterval   = configuration->plotInterval;

    if (Dimensions == 3) {
      _simulation->_memory.parameters()->g(2) =
        configuration->environment(2);
    }

    _simulation->_memory.initialize(rank,
                                    processor_size,
                                    global_cell_size,
                                    geometry_width);

    setLeftWallAs(configuration->walls[0][0]->type());
    setRightWallAs(configuration->walls[0][1]->type());
    setBottomWallAs(configuration->walls[1][0]->type());
    setTopWallAs(configuration->walls[1][1]->type());
    setBackWallAs(configuration->walls[2][0]->type());
    setFrontWallAs(configuration->walls[2][1]->type());
  }

  SolverType*
  simulation() const {
    return _simulation;
  }

  void
  setLeftWallAs(FluidSimulation::WallEnum type) {
    GhostHandlersBuilderType<0, 0> handlers(_configuration, _simulation);

    if (_simulation->_memory.parallelDistribution()->neighbors[0][0] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == FluidSimulation::WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == FluidSimulation::WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == FluidSimulation::WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setRightWallAs(FluidSimulation::WallEnum type) {
    GhostHandlersBuilderType<0, 1> handlers(_configuration, _simulation);

    if (_simulation->_memory.parallelDistribution()->neighbors[0][1] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == FluidSimulation::WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == FluidSimulation::WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == FluidSimulation::WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setBottomWallAs(FluidSimulation::WallEnum type) {
    if (Dimensions < 2) { return; }

    GhostHandlersBuilderType<1, 0> handlers(_configuration, _simulation);

    if (_simulation->_memory.parallelDistribution()->neighbors[1][0] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == FluidSimulation::WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == FluidSimulation::WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == FluidSimulation::WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setTopWallAs(FluidSimulation::WallEnum type) {
    if (Dimensions < 2) { return; }

    GhostHandlersBuilderType<1, 1> handlers(_configuration, _simulation);

    if (_simulation->_memory.parallelDistribution()->neighbors[1][1] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == FluidSimulation::WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == FluidSimulation::WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == FluidSimulation::WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setBackWallAs(FluidSimulation::WallEnum type) {
    if (Dimensions < 3) { return; }

    GhostHandlersBuilderType<2, 0> handlers(_configuration, _simulation);

    if (_simulation->_memory.parallelDistribution()->neighbors[2][0] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == FluidSimulation::WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == FluidSimulation::WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == FluidSimulation::WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

  void
  setFrontWallAs(FluidSimulation::WallEnum type) {
    if (Dimensions < 3) { return; }

    GhostHandlersBuilderType<2, 1> handlers(_configuration, _simulation);

    if (_simulation->_memory.parallelDistribution()->neighbors[2][1] >= 0) {
      handlers.setAsMpiExchange();
    } else if (type == FluidSimulation::WallEnum::Input) {
      handlers.setAsInput();
    } else if (type == FluidSimulation::WallEnum::ParabolicInput) {
      handlers.setAsParabolicInput();
    } else if (type == FluidSimulation::WallEnum::Output) {
      handlers.setAsOutput();
    }
  }

private:
  FluidSimulation::Configuration* _configuration;
  SolverType*                     _simulation;
};
}
}
