#pragma once

#include "GhostLayer/IfsfdHandlers.hpp"
#include "GhostLayer/SfsfdHandlers.hpp"
#include "ImmersedBoundary/Controller.hpp"
#include "PeSolver.hpp"
#include "SolverTraits.hpp"
#include "IfsfdMemory.hpp"
#include "SfsfdMemory.hpp"

#include <Uni/Logging/macros>

namespace precice {
class SolverInterface;
}

namespace FsiSimulation {
namespace FluidSimulation {
class Reporter;
template <typename TSolverTraits>
class FsfdSolver {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions,
    SolverId   = SolverTraitsType::SolverId
  };

  using GridGeometryType = typename SolverTraitsType::GridGeometryType;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  using ImmersedBoundaryControllerType
          = typename SolverTraitsType::ImmersedBoundaryControllerType;

  using GhostHandlersType
          = typename SolverTraitsType::GhostHandlersType;

  using PeSolverType
          = typename SolverTraitsType::PeSolverType;

public:
  FsfdSolver() {}

  FsfdSolver(FsfdSolver const& other) = delete;

  ~FsfdSolver() {}

  FsfdSolver const&
  operator=(FsfdSolver const& other) = delete;

  MemoryType const*
  memory() const {
    return &_memory;
  }

  MemoryType*
  memory() {
    return &_memory;
  }

  ImmersedBoundaryControllerType const*
  immersedBoundaryController() const {
    return &_ibController;
  }

  ImmersedBoundaryControllerType*
  immersedBoundaryController() {
    return &_ibController;
  }

  GhostHandlersType const*
  ghostHandlers() const {
    return &_ghostHandlers;
  }

  GhostHandlersType*
  ghostHandlers() {
    return &_ghostHandlers;
  }

  void
  initialize(precice::SolverInterface*, Reporter*);

  void
  iterate();

private:
  MemoryType                     _memory;
  PeSolverType                   _peSolver;
  GhostHandlersType              _ghostHandlers;
  ImmersedBoundaryControllerType _ibController;

  Reporter* _reporter;
};

extern template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 0, double, 2 >>;
extern template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 1, double, 2 >>;
extern template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 0, double, 2 >>;
extern template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 1, double, 2 >>;
extern template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 0, double, 3 >>;
extern template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 1, double, 3 >>;
extern template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 0, double, 3 >>;
extern template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 1, double, 3 >>;

extern template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 0, double, 2 >>;
extern template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 1, double, 2 >>;
extern template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 0, double, 2 >>;
extern template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 1, double, 2 >>;
extern template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 0, double, 3 >>;
extern template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 1, double, 3 >>;
extern template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 0, double, 3 >>;
extern template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 1, double, 3 >>;
}
}
