#pragma once

#include "FsfdSolver.hpp"
#include "Grid.hpp"
#include "IfsfdCellAccessor.hpp"
#include "IfsfdMemory.hpp"
#include "ParallelDistribution.hpp"
#include "Parameters.hpp"
#include "SfsfdCellAccessor.hpp"
#include "SfsfdMemory.hpp"
#include "PeSolver.hpp"

#include "GhostLayer/IfsfdHandlers.hpp"
#include "GhostLayer/SfsfdHandlers.hpp"

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TGridGeometry,
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

  using MemoryType = SfsfdMemory<SfsfdSolverTraits>;

  using CellAccessorType = SfsfdCellAccessor<SfsfdSolverTraits>;

  using BaseGridType = Uni::StructuredGrid::Basic::Grid<CellAccessorType>;

  using ParallelDistributionType = ParallelDistribution<Dimensions>;

  using VectorDsType = Eigen::Matrix<ScalarType, Dimensions, 1>;

  using VectorDiType = Eigen::Matrix<int, Dimensions, 1>;

  using GhostHandlersType = typename GhostLayer::SfsfdHandlers<Dimensions>;

  using PeSolverType = SfsfdPeSolver<SfsfdSolverTraits>;

  using SolverType = FsfdSolver<SfsfdSolverTraits, 0>;
};

template <typename TGridGeometry,
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

  using MemoryType = IfsfdMemory<IfsfdSolverTraits>;

  using CellAccessorType = IfsfdCellAccessor<IfsfdSolverTraits>;

  using BaseGridType = Uni::StructuredGrid::Basic::Grid<CellAccessorType>;

  using ParallelDistributionType = ParallelDistribution<Dimensions>;

  using VectorDsType = Eigen::Matrix<ScalarType, Dimensions, 1>;

  using VectorDiType = Eigen::Matrix<int, Dimensions, 1>;

  using GhostHandlersType = typename GhostLayer::IfsfdHandlers<Dimensions>;

  using PeSolverType = IfsfdPeSolver<IfsfdSolverTraits>;

  using SolverType = FsfdSolver<IfsfdSolverTraits, 1>;
};
}
}
