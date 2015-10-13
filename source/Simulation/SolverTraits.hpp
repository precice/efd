#pragma once

#include <Eigen/Core>

#include <Uni/StructuredGrid/Basic/Grid>

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
class IfsfdMemory;
template <typename T, unsigned TSolverId = T::SolverId>
class PeSolver;
template <typename T>
class FsfdSolver;
namespace GhostLayer {
template <int D>
class SfsfdHandlers;
template <int D>
class IfsfdHandlers;
}

template <unsigned TDebug,
          typename TScalar,
          unsigned TDimensions>
struct SfsfdSolverTraits {
  using ScalarType = TScalar;

  enum {
    Dimensions = TDimensions,
    SolverId   = 0,
    DebugLevel = TDebug
  };

  using GridType = Grid<SfsfdSolverTraits>;

  using ParametersType = Parameters<ScalarType, Dimensions>;

  using GridGeometryType = UniformGridGeometry<ScalarType, Dimensions>;

  using MemoryType = SfsfdMemory<SfsfdSolverTraits>;

  using CellAccessorType = SfsfdCellAccessor<SfsfdSolverTraits>;

  using SubgridType = Uni::StructuredGrid::Basic::Grid<
          Uni::StructuredGrid::Basic::DummyForThisGrid,
          CellAccessorType,
          Dimensions>;

  using ParallelDistributionType = ParallelDistribution<Dimensions>;

  using VectorDsType = Eigen::Matrix<ScalarType, Dimensions, 1>;

  using VectorDiType = Eigen::Matrix<int, Dimensions, 1>;

  using GhostHandlersType = typename GhostLayer::SfsfdHandlers<Dimensions>;

  using PeSolverType = PeSolver<SfsfdSolverTraits>;

  using SolverType = FsfdSolver<SfsfdSolverTraits>;
};

template <int TDebug,
          typename TScalar,
          unsigned TDimensions>
struct IfsfdSolverTraits {
  using ScalarType = TScalar;

  enum {
    Dimensions = TDimensions,
    SolverId   = 1,
    DebugLevel = TDebug
  };

  using GridType = Grid<IfsfdSolverTraits>;

  using ParametersType = Parameters<ScalarType, Dimensions>;

  using GridGeometryType = UniformGridGeometry<ScalarType, Dimensions>;

  using MemoryType = IfsfdMemory<IfsfdSolverTraits>;

  using CellAccessorType = IfsfdCellAccessor<IfsfdSolverTraits>;

  using SubgridType = Uni::StructuredGrid::Basic::Grid<
          Uni::StructuredGrid::Basic::DummyForThisGrid,
          CellAccessorType,
          Dimensions>;

  using ParallelDistributionType = ParallelDistribution<Dimensions>;

  using VectorDsType = Eigen::Matrix<ScalarType, Dimensions, 1>;

  using VectorDiType = Eigen::Matrix<int, Dimensions, 1>;

  using GhostHandlersType = typename GhostLayer::IfsfdHandlers<Dimensions>;

  using PeSolverType = PeSolver<IfsfdSolverTraits>;

  using SolverType = FsfdSolver<IfsfdSolverTraits>;
};
}
}

#define Fluid_DeclareExternTemplates(ClassName) \
  extern template class ClassName               \
    < SfsfdSolverTraits < 0, double, 2 >>;      \
  extern template class ClassName               \
    < SfsfdSolverTraits < 1, double, 2 >>;      \
  extern template class ClassName               \
    < SfsfdSolverTraits < 0, double, 3 >>;      \
  extern template class ClassName               \
    < SfsfdSolverTraits < 1, double, 3 >>;      \
  extern template class ClassName               \
    < IfsfdSolverTraits < 0, double, 2 >>;      \
  extern template class ClassName               \
    < IfsfdSolverTraits < 1, double, 2 >>;      \
  extern template class ClassName               \
    < IfsfdSolverTraits < 0, double, 3 >>;      \
  extern template class ClassName               \
    < IfsfdSolverTraits < 1, double, 3 >>;

#define Fluid_InstantiateExternTemplates(ClassName) \
  template class ClassName                          \
    < SfsfdSolverTraits < 0, double, 2 >>;          \
  template class ClassName                          \
    < SfsfdSolverTraits < 1, double, 2 >>;          \
  template class ClassName                          \
    < SfsfdSolverTraits < 0, double, 3 >>;          \
  template class ClassName                          \
    < SfsfdSolverTraits < 1, double, 3 >>;          \
  template class ClassName                          \
    < IfsfdSolverTraits < 0, double, 2 >>;          \
  template class ClassName                          \
    < IfsfdSolverTraits < 1, double, 2 >>;          \
  template class ClassName                          \
    < IfsfdSolverTraits < 0, double, 3 >>;          \
  template class ClassName                          \
    < IfsfdSolverTraits < 1, double, 3 >>;
