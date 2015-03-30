#pragma once

#include "ImmersedBoundary/Controller.hpp"
#include "LinearSolver.hpp"
#include "Private/mpigenerics.hpp"
#include "functions.hpp"

#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits, int TSolverType>
class FsfdSolver {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
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
  initialize(precice::SolverInterface* preciceInteface) {
    _ibController.initialize(preciceInteface);

    _peSolver.initialize(&_memory, &_ghostHandlers);

    for (auto& accessor : _memory.grid()->innerGrid) {
      _ibController.computePositionInRespectToGeometry(accessor);
      accessor.velocity() = VectorDsType::Zero();
      accessor.setDiffusion(VectorDsType::Zero());
      accessor.setForce(VectorDsType::Zero());
      accessor.setBodyForce(VectorDsType::Zero());
      PressureProcessing<TSolverType>::initialize(accessor, 0.0);
      computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
        (accessor, _memory.maxVelocity());
    }

    _ghostHandlers.executeVelocityInitialization();
    _ghostHandlers.executeVelocityMpiExchange();
    _ghostHandlers.executePressureMpiExchange();

    for (auto const& accessor : _memory.grid()->innerGrid) {
      _ibController.createFluidMeshVertex(accessor);

      auto convection
        = ConvectionProcessing<Dimensions>::compute(accessor,
                                                    _memory.parameters()->gamma());
      accessor.convection() = convection;
    }

    _memory.maxVelocity()     = VectorDsType::Zero();
    _memory.timeStepSize()    = 1.0;
    _memory.time()            = 0.0;
    _memory.iterationNumber() = 0;
  }

  void
  iterate() {
    _memory.timeStepSize()
      = TimeStepProcessing<ScalarType, Dimensions>::compute(
      _memory.parameters()->re(),
      _memory.parameters()->tau(),
      _memory.gridGeometry()->minCellWidth(),
      _memory.maxVelocity());

    _memory.maxVelocity() = VectorDsType::Zero();

    for (auto const& accessor : _memory.grid()->innerGrid) {
      _memory.setForceAt(accessor.globalIndex(), VectorDsType::Zero());

      auto convection = ConvectionProcessing<Dimensions>::compute(
        accessor,
        _memory.parameters()->gamma());

      auto previousConvection = accessor.convection();

      if (_memory.iterationNumber() == 0) {
        previousConvection = convection;
      }

      accessor.convection() = convection;

      auto diffusion = DiffusionProcessing<Dimensions>::compute(accessor);

      accessor.setDiffusion(diffusion);

      diffusion = diffusion / _memory.parameters()->re();

      VectorDsType velocity;

      for (int d = 0; d < Dimensions; ++d) {
        velocity(d) = 0.5 * (accessor.velocity(d, -1, d) +
                             accessor.velocity(d));
      }

      VectorDsType grad_pressure;

      for (int d = 0; d < Dimensions; ++d) {
        grad_pressure(d)
          = (accessor.pressure(d, +1) - accessor.pressure())
            / (0.5 * (accessor.width(d, +1, d) + accessor.width(d)));
      }

      if (TSolverType == 0) {
        accessor.fgh()
          = accessor.velocity()
            + _memory.timeStepSize() * (diffusion
                                        - convection
                                        + _memory.parameters()->g()
                                        );
      } else {
        accessor.fgh()
          = accessor.velocity()
            + _memory.timeStepSize() * (-1.5 * convection
                                        + 0.5 * previousConvection
                                        + 0.5 * diffusion
                                        - grad_pressure
                                        + _memory.parameters()->g());
      }

      auto status = _ibController.doesVertexExist(accessor);

      if (status.second) {
        VectorDsType temp
          = velocity
            + _memory.timeStepSize() * (-1.5 * convection
                                        + 0.5 * previousConvection
                                        + diffusion
                                        - grad_pressure);

        _ibController.writeFluidVelocity(status.first, temp);
      }
    }

    _ibController.mapData(_memory.timeStepSize());

    VectorDsType total_force = VectorDsType::Zero();

    for (auto it = _ibController.begin(); it != _ibController.end(); ++it) {
      VectorDsType force;
      _ibController.readFluidForce(it, force);

      _memory.fgh()[it->first] += force;
      _memory.setForceAt(it->first, force);
      total_force += force;
    }
    logInfo("{1}", total_force.transpose());

    _peSolver.executeVpe();
    _ghostHandlers.executeFghMpiExchange();

    _peSolver.executePpe();
    _ghostHandlers.executePressureMpiExchange();

    for (auto& accessor : _memory.grid()->innerGrid) {
      VectorDsType grad_pressure
        = PressureProcessing<TSolverType>::grad(accessor);

      accessor.velocity()
        = accessor.fgh()
          - _memory.timeStepSize() * grad_pressure;

      computeMaxVelocity<CellAccessorType, ScalarType, Dimensions>
        (accessor, _memory.maxVelocity());
    }

    _ghostHandlers.executeVelocityInitialization();
    _ghostHandlers.executeVelocityMpiExchange();

    _ibController.computeBodyForce(&_memory);

    _memory.time() += _memory.timeStepSize();
    ++_memory.iterationNumber();
  }

private:
  MemoryType                     _memory;
  PeSolverType                   _peSolver;
  GhostHandlersType              _ghostHandlers;
  ImmersedBoundaryControllerType _ibController;
};
}
}
