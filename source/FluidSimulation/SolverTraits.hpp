#pragma once

#include "FsfdSolver.hpp"
#include "Grid.hpp"
#include "IfsfdCellAccessor.hpp"
#include "IfsfdMemory.hpp"
#include "ImmersedBoundary/Controller.hpp"
#include "ParallelDistribution.hpp"
#include "Parameters.hpp"
#include "PeSolver.hpp"
#include "SfsfdCellAccessor.hpp"
#include "SfsfdMemory.hpp"

#include "GhostLayer/IfsfdHandlers.hpp"
#include "GhostLayer/SfsfdHandlers.hpp"

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits,
          int TSolverType,
          int TDebug>
struct MemoryTraits {};

template <typename TSolverTraits>
struct MemoryTraits<TSolverTraits, 0, 0> {
  using Type = SfsfdMemory<TSolverTraits>;
};

template <typename TSolverTraits>
struct MemoryTraits<TSolverTraits, 0, 1> {
  using Type = SfsfdDebugMemory<TSolverTraits>;
};

template <typename TSolverTraits>
struct MemoryTraits<TSolverTraits, 1, 0> {
  using Type = IfsfdMemory<TSolverTraits>;
};

template <typename TSolverTraits>
struct MemoryTraits<TSolverTraits, 1, 1> {
  using Type = IfsfdDebugMemory<TSolverTraits>;
};

template <typename TSolverTraits,
          int TImmersedBoudnaryType>
struct ImmersedBoundaryControllerTraits {};

template <typename TSolverTraits>
struct ImmersedBoundaryControllerTraits<TSolverTraits, 0> {
  using Type = ImmersedBoundary::BasicController<TSolverTraits>;
};

template <typename TSolverTraits>
struct ImmersedBoundaryControllerTraits<TSolverTraits, 1> {
  using Type = ImmersedBoundary::Controller<TSolverTraits>;
};

template <typename TGridGeometry,
          int TImmersedBoudnaryType,
          int TDebug,
          typename TScalar,
          int TD>
struct SfsfdSolverTraits {
  using ScalarType = TScalar;

  enum {
    Dimensions = TD
  };

  using GridType = Grid<SfsfdSolverTraits>;

  using ParametersType = Parameters<ScalarType, Dimensions>;

  using GridGeometryType = TGridGeometry;

  using MemoryType
          = typename MemoryTraits<SfsfdSolverTraits, 0, TDebug>::Type;

  using CellAccessorType = SfsfdCellAccessor<SfsfdSolverTraits>;

  using BaseGridType = Uni::StructuredGrid::Basic::Grid<CellAccessorType>;

  using ParallelDistributionType = ParallelDistribution<Dimensions>;

  using VectorDsType = Eigen::Matrix<ScalarType, Dimensions, 1>;

  using VectorDiType = Eigen::Matrix<int, Dimensions, 1>;

  using GhostHandlersType = typename GhostLayer::SfsfdHandlers<Dimensions>;

  using PeSolverType = SfsfdPeSolver<SfsfdSolverTraits>;

  using ImmersedBoundaryControllerType
          = typename ImmersedBoundaryControllerTraits
            <SfsfdSolverTraits, TImmersedBoudnaryType>::Type;

  using SolverType = FsfdSolver<SfsfdSolverTraits, 0>;
};

template <typename TGridGeometry,
          int TImmersedBoudnaryType,
          int TDebug,
          typename TScalar,
          int TD>
struct IfsfdSolverTraits {
  using ScalarType = TScalar;

  enum {
    Dimensions = TD
  };

  using GridType = Grid<IfsfdSolverTraits>;

  using ParametersType = Parameters<ScalarType, Dimensions>;

  using GridGeometryType = TGridGeometry;

  using MemoryType
          = typename MemoryTraits<IfsfdSolverTraits, 1, TDebug>::Type;

  using CellAccessorType = IfsfdCellAccessor<IfsfdSolverTraits>;

  using BaseGridType = Uni::StructuredGrid::Basic::Grid<CellAccessorType>;

  using ParallelDistributionType = ParallelDistribution<Dimensions>;

  using VectorDsType = Eigen::Matrix<ScalarType, Dimensions, 1>;

  using VectorDiType = Eigen::Matrix<int, Dimensions, 1>;

  using GhostHandlersType = typename GhostLayer::IfsfdHandlers<Dimensions>;

  using PeSolverType = IfsfdPeSolver<IfsfdSolverTraits>;

  using ImmersedBoundaryControllerType
          = typename ImmersedBoundaryControllerTraits
            <IfsfdSolverTraits, TImmersedBoudnaryType>::Type;

  using SolverType = FsfdSolver<IfsfdSolverTraits, 1>;
};
}
}
