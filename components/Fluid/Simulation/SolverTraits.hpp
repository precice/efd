#pragma once

#include <Eigen/Core>

namespace Uni {
namespace StructuredGrid {
namespace Basic {
template <typename T>
class Grid;
}
}
}
namespace FsiSimulation {
namespace FluidSimulation {
template <unsigned D>
class ParallelDistribution;
template <typename T, unsigned D>
class Parameters;
template <typename T, unsigned D>
class UniformGridGeometry;
template <typename T>
class Grid;
template <typename T>
class SfsfdCellAccessor;
template <typename T>
class IfsfdCellAccessor;
template <typename T>
class SfsfdMemory;
template <typename T>
class SfsfdDebugMemory;
template <typename T>
class IfsfdMemory;
template <typename T>
class IfsfdDebugMemory;
template <typename T>
class SfsfdPeSolver;
template <typename T>
class IfsfdPeSolver;
template <typename T>
class FsfdSolver;
namespace ImmersedBoundary {
template <typename T>
class BasicController;
template <typename T>
class Controller;
}
namespace GhostLayer {
template <int D>
class SfsfdHandlers;
template <int D>
class IfsfdHandlers;
}

template <typename TSolverTraits,
          int TSolverType,
          unsigned TDebug>
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

template <int TImmersedBoudnaryType,
          unsigned TDebug,
          typename TScalar,
          unsigned TDimensions>
struct SfsfdSolverTraits {
  using ScalarType = TScalar;

  enum {
    Dimensions = TDimensions,
    SolverId   = 0
  };

  using GridType = Grid<SfsfdSolverTraits>;

  using ParametersType = Parameters<ScalarType, Dimensions>;

  using GridGeometryType = UniformGridGeometry<ScalarType, Dimensions>;

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

  using SolverType = FsfdSolver<SfsfdSolverTraits>;
};

template <int TImmersedBoudnaryType,
          int TDebug,
          typename TScalar,
          unsigned TDimensions>
struct IfsfdSolverTraits {
  using ScalarType = TScalar;

  enum {
    Dimensions = TDimensions,
    SolverId   = 1
  };

  using GridType = Grid<IfsfdSolverTraits>;

  using ParametersType = Parameters<ScalarType, Dimensions>;

  using GridGeometryType = UniformGridGeometry<ScalarType, Dimensions>;

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

  using SolverType = FsfdSolver<IfsfdSolverTraits>;
};
}
}
