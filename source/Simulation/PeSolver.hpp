#pragma once

#include "LinearSolver.hpp"
#include "SolverTraits.hpp"

#include <memory>

namespace Fluid {
namespace Simulation {
template <typename TSolverTraits,
          unsigned TSolverId>
class PeSolver {};

template <typename TSolverTraits>
class PeSolver<TSolverTraits, 0> {
public:
  using SolverTraitsType = TSolverTraits;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using GhostHandlersType = typename SolverTraitsType::GhostHandlersType;

  void
  initialize(MemoryType const*        memory,
             GhostHandlersType const* ghost_handlers);

  void
  executeVpe();

  void
  executePpe();

private:
  GhostHandlersType const*      _ghostHandlers;
  std::unique_ptr<LinearSolver> _ppeSolver;
};

namespace Private {
template <typename TSolverTraits, int TDimensions>
class IfsfdPeSolver3DPart {};

template <typename TSolverTraits>
class IfsfdPeSolver3DPart<TSolverTraits, 2> {
public:
  using SolverTraitsType = TSolverTraits;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using GhostHandlersType = typename SolverTraitsType::GhostHandlersType;

  void
  initialize(MemoryType const*, GhostHandlersType const*) {}

  void
  executeVpe() {}
};

template <typename TSolverTraits>
class IfsfdPeSolver3DPart<TSolverTraits, 3> {
public:
  using SolverTraitsType = TSolverTraits;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using GhostHandlersType = typename SolverTraitsType::GhostHandlersType;

  void
  initialize(MemoryType const*, GhostHandlersType const*);

  void
  executeVpe();

private:
  MemoryType const*             _memory;
  long double                   _vzpeAbsoluteTolerance;
  long double                   _vzpeRelativeTolerance;
  std::unique_ptr<LinearSolver> _vzpeSolver;
};
}

template <typename TSolverTraits>
class PeSolver<TSolverTraits, 1> {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using GhostHandlersType = typename SolverTraitsType::GhostHandlersType;

  using VzpeSolverType
          = Private::IfsfdPeSolver3DPart<SolverTraitsType, Dimensions>;

  void
  initialize(MemoryType const*        memory,
             GhostHandlersType const* ghost_handlers);

  void
  executeVpe();

  void
  executePpe();

private:
  MemoryType const* _memory;

  long double                   _vxpeAbsoluteTolerance;
  long double                   _vxpeRelativeTolerance;
  std::unique_ptr<LinearSolver> _vxpeSolver;

  long double                   _vypeAbsoluteTolerance;
  long double                   _vypeRelativeTolerance;
  std::unique_ptr<LinearSolver> _vypeSolver;

  VzpeSolverType _vzpeSolver;

  std::unique_ptr<LinearSolver> _ppeSolver;
};

Fluid_DeclareExternTemplates(PeSolver);
}
}
