#pragma once

#include "SolverTraits.hpp"

#include <Uni/Firewall/Implementation>

namespace precice {
class SolverInterface;
}

namespace FsiSimulation {
namespace FluidSimulation {
class Configuration;

namespace ImmersedBoundary {
template <typename TVector>
class Controller;
}

template <typename T>
class FsfdSolverImplementation;
class Reporter;
template <typename TSolverTraits>
class FsfdSolver {
private:
  using Implementation = FsfdSolverImplementation<TSolverTraits>;

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

  using IbControllerType
          = ImmersedBoundary::Controller<VectorDsType>;

  using GhostHandlersType
          = typename SolverTraitsType::GhostHandlersType;

  using PeSolverType
          = typename SolverTraitsType::PeSolverType;

public:
  FsfdSolver(Configuration const* configuration);

  FsfdSolver(FsfdSolver const& other) = delete;

  ~FsfdSolver();

  FsfdSolver const&
  operator=(FsfdSolver const& other) = delete;

  MemoryType const*
  memory() const;

  MemoryType*
  memory();

  GhostHandlersType const*
  ghostHandlers() const;

  GhostHandlersType*
  ghostHandlers();

  void
  initialize(precice::SolverInterface* precice_interface,
             Reporter*                 reporter);

  ScalarType
  computeTimeStepSize();

  void
  iterate();

private:
  void
  iterateWithFastIbVelocityPrediction();

  void
  iterateWithFullIbVelocityPrediction();

  std::array<VectorDsType, 3>
  updateFgh(CellAccessorType const& accessor);

  void
  advanceFsi();

  void
  addIbForces();

  void
  solvePoissonEquations();

  void
  updateVelocities();

  void
  finalizeIteration();

  void
  locateStructure();

  void
  sendCouplingData();

  void
  computeBodyForce();

  Uni_Firewall_IMPLEMENTATION_LINK(FsfdSolverImplementation<TSolverTraits> );
};

Fluid_DeclareExternTemplates(FsfdSolver);
}
}
