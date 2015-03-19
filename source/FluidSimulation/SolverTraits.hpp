#pragma once

#include "CellAccessor.hpp"
#include "Grid.hpp"
#include "Memory.hpp"
#include "ParallelDistribution.hpp"
#include "Parameters.hpp"

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TGridGeometry,
          typename TScalar,
          int TD>
struct SolverTraits {
  using ScalarType = TScalar;

  enum {
    Dimensions = TD
  };

  using GridType = Grid<SolverTraits>;

  using ParametersType = Parameters<ScalarType, Dimensions>;

  using GridGeometryType = TGridGeometry;

  using MemoryType = Memory<SolverTraits>;

  using CellAccessorType = CellAccessor<SolverTraits>;

  using BaseGridType = Uni::StructuredGrid::Basic::Grid<CellAccessorType>;

  using ParallelDistributionType = ParallelDistribution<Dimensions>;

  using VectorDsType = Eigen::Matrix<ScalarType, Dimensions, 1>;

  using VectorDiType = Eigen::Matrix<int, Dimensions, 1>;
};
}
}
