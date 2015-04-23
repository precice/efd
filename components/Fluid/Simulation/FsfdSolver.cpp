#include "FsfdSolver.hpp"

#include "GhostLayer/IfsfdHandlers.hpp"
#include "GhostLayer/SfsfdHandlers.hpp"
#include "Grid.hpp"
#include "GridGeometry.hpp"
#include "IfsfdCellAccessor.hpp"
#include "IfsfdMemory.hpp"
#include "ImmersedBoundary/Controller.hpp"
#include "PeSolver.hpp"
#include "Private/mpigenerics.hpp"
#include "Reporter.hpp"
#include "SfsfdCellAccessor.hpp"
#include "SfsfdMemory.hpp"
#include "computemaxvelocity.hpp"
#include "functions.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/Firewall/Interface>
#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TScalar, typename TVector>
inline TScalar
compute_time_step_size(TScalar const& re,
                       TScalar const& tau,
                       TVector const& minCellWidth,
                       TVector const& maxVelocity) {
  TScalar factor = minCellWidth.cwiseProduct(minCellWidth).cwiseInverse().sum();

  TScalar localMin = (TScalar)(re / (2.0 * factor));

  for (unsigned d = 0; d < TVector::RowsAtCompileTime; ++d) {
    if (std::abs(maxVelocity(d)) > std::numeric_limits<TScalar>::epsilon()) {
      localMin = std::min(localMin,
                          1.0 / maxVelocity(d));
    }
  }

  TScalar globalMin = std::numeric_limits<TScalar>::max();
  Private::mpiAllReduce<TScalar>(&localMin,
                                 &globalMin,
                                 1,
                                 MPI_MIN,
                                 PETSC_COMM_WORLD);

  factor  = globalMin;
  factor *= tau;

  return factor;
}

template <typename TSolverTraits>
class FsfdSolverImplementation {
public:
  using Interface = FsfdSolver<TSolverTraits>;

  FsfdSolverImplementation(Interface* in) : _in(in) {}

  typename Interface::MemoryType                     memory;
  typename Interface::PeSolverType                   peSolver;
  typename Interface::GhostHandlersType              ghostHandlers;
  typename Interface::ImmersedBoundaryControllerType ibController;
  Reporter* reporter;

  Uni_Firewall_INTERFACE_LINK(FsfdSolver<TSolverTraits> );
};

template <typename T>
FsfdSolver<T>::
FsfdSolver() : _im(new Implementation(this)) {}

template <typename T>
FsfdSolver<T>::
~FsfdSolver() {}

template <typename T>
typename FsfdSolver<T>::MemoryType const*
FsfdSolver<T>::
memory() const {
  return &_im->memory;
}

template <typename T>
typename FsfdSolver<T>::MemoryType *
FsfdSolver<T>::
memory() {
  return &_im->memory;
}

template <typename T>
typename FsfdSolver<T>::ImmersedBoundaryControllerType const*
FsfdSolver<T>::
immersedBoundaryController() const {
  return &_im->ibController;
}

template <typename T>
typename FsfdSolver<T>::ImmersedBoundaryControllerType *
FsfdSolver<T>::
immersedBoundaryController() {
  return &_im->ibController;
}

template <typename T>
typename FsfdSolver<T>::GhostHandlersType const*
FsfdSolver<T>::
ghostHandlers() const {
  return &_im->ghostHandlers;
}

template <typename T>
typename FsfdSolver<T>::GhostHandlersType *
FsfdSolver<T>::
ghostHandlers() {
  return &_im->ghostHandlers;
}

template <typename T>
void
FsfdSolver<T>::
initialize(precice::SolverInterface* preciceInteface,
           Reporter*                 reporter) {
  _im->ibController.initialize(preciceInteface, &_im->memory);
  _im->reporter = reporter;

  _im->peSolver.initialize(&_im->memory, &_im->ghostHandlers);

  _im->memory.maxVelocity()     = VectorDsType::Zero();
  _im->memory.timeStepSize()    = 0.0;
  _im->memory.time()            = 0.0;
  _im->memory.iterationNumber() = 0;

  for (auto& accessor : _im->memory.grid()->innerGrid) {
    _im->ibController.computePositionInRespectToGeometry(accessor);
    accessor.velocity() = VectorDsType::Zero();
    accessor.setDiffusion(VectorDsType::Zero());
    accessor.setForce(VectorDsType::Zero());
    accessor.setBodyForce(VectorDsType::Zero());
    PressureProcessing<SolverId>::initialize(accessor, 0.0);
    compute_max_velocity(accessor, _im->memory.maxVelocity());
  }

  _im->ghostHandlers.executeVelocityInitialization();
  _im->ghostHandlers.executeVelocityMpiExchange();
  _im->ghostHandlers.executePressureMpiExchange();

  for (auto const& accessor : _im->memory.grid()->innerGrid) {
    _im->ibController.createFluidMeshVertex(accessor);

    auto convection
      = ConvectionProcessing<Dimensions>::compute(accessor,
                                                  _im->memory.parameters()->
                                                  gamma());
    accessor.convection() = convection;
  }

  _im->reporter->setAt("Re",
                       _im->memory.parameters()->re());
  _im->reporter->setAt("Gamma",
                       _im->memory.parameters()->gamma());
  _im->reporter->setAt("Tau",
                       _im->memory.parameters()->tau());
  _im->reporter->setAt("G",
                       _im->memory.parameters()->g());
  _im->reporter->setAt("OuterLayerSize",
                       _im->ibController.outerLayerSize());
  _im->reporter->setAt("InnerLayerSize",
                       _im->ibController.innerLayerSize());
  _im->reporter->setAt("ProcessorSize",
                       _im->memory.parallelDistribution()->processorSize);
  _im->reporter->setAt("GlobalCellSize",
                       _im->memory.parallelDistribution()->globalCellSize);
  _im->reporter->setAt("UniformLocalCellSize",
                       _im->memory.parallelDistribution()->uniformLocalCellSize);
  _im->reporter->setAt("LastLocalCellSize",
                       _im->memory.parallelDistribution()->lastLocalCellSize);
  _im->reporter->setAt("Width",
                       _im->memory.gridGeometry()->size());
  _im->reporter->setAt("CellWidth",
                       _im->memory.gridGeometry()->minCellWidth());
  _im->reporter->setAt("SolverId", SolverId);

  _im->reporter->addAt(0, "IterationNumber");
  _im->reporter->addAt(1, "Time");
  _im->reporter->addAt(2, "TimeStepSize");
  _im->reporter->addAt(3, "Force1");
  _im->reporter->addAt(4, "Force2");
  _im->reporter->addAt(5, "Force3");

  _im->reporter->recordIteration();

  _im->reporter->addAt(0, _im->memory.iterationNumber());
  _im->reporter->addAt(1, _im->memory.time());
  _im->reporter->addAt(2, _im->memory.timeStepSize());
  _im->reporter->addAt(3, VectorDsType::Zero().eval());
  _im->reporter->addAt(4, VectorDsType::Zero().eval());
  _im->reporter->addAt(5, VectorDsType::Zero().eval());

  _im->reporter->recordIteration();
}

template <typename T>
void
FsfdSolver<T>::
iterate() {
  _im->memory.timeStepSize()
    = compute_time_step_size(_im->memory.parameters()->re(),
                             _im->memory.parameters()->tau(),
                             _im->memory.gridGeometry()->minCellWidth(),
                             _im->memory.maxVelocity());

  if (_im->memory.parallelDistribution()->rank == 0) {
    logInfo("dt = {1}", _im->memory.timeStepSize());
    logInfo("maxv = {1}",
            _im->memory.maxVelocity().cwiseProduct(
              _im->memory.gridGeometry()->minCellWidth()).transpose());
  }
  _im->memory.maxVelocity()
    = VectorDsType::Constant(std::numeric_limits<ScalarType>::min());

  _im->ibController.resetIntermediateData();

  for (auto const& accessor : _im->memory.grid()->innerGrid) {
    _im->memory.setForceAt(accessor.globalIndex(), VectorDsType::Zero());

    auto convection = ConvectionProcessing<Dimensions>::compute(
      accessor,
      _im->memory.parameters()->gamma());

    auto previousConvection = accessor.convection();

    if (_im->memory.iterationNumber() == 0) {
      previousConvection = convection;
    }

    accessor.convection() = convection;

    auto diffusion = DiffusionProcessing<Dimensions>::compute(accessor);

    accessor.setDiffusion(diffusion);

    diffusion = diffusion / _im->memory.parameters()->re();

    VectorDsType velocity;

    // for (int d = 0; d < Dimensions; ++d) {
    // velocity(d) = 0.5 * (accessor.velocity(d, -1, d) +
    // accessor.velocity(d));
    // }

    VectorDsType grad_pressure;

    for (int d = 0; d < Dimensions; ++d) {
      grad_pressure(d)
        = (accessor.pressure(d, +1) - accessor.pressure())
          / (0.5 * (accessor.width(d, +1, d) + accessor.width(d)));
    }

    if (SolverId == 0) {
      accessor.fgh()
        = accessor.velocity()
          + _im->memory.timeStepSize() * (diffusion
                                          - convection
                                          + _im->memory.parameters()->g()
                                          );
    } else {
      accessor.fgh()
        = accessor.velocity()
          + _im->memory.timeStepSize() * (-1.5 * convection
                                          + 0.5 * previousConvection
                                          + 0.5 * diffusion
                                          - grad_pressure
                                          + _im->memory.parameters()->g());
      // logInfo("Fgh | {1}", accessor.index().transpose());
    }

    auto status = _im->ibController.doesVertexExist(accessor);

    if (status.second) {
      _im->memory.setForceAt(accessor.globalIndex(), VectorDsType::Zero());
      auto temp
        = accessor.velocity()
          + _im->memory.timeStepSize() * (-1.5 * convection
                                          + 0.5 * previousConvection
                                          + diffusion
                                          - grad_pressure);

      _im->ibController.setFluidVelocity(status.first, temp);
    }
  }

  _im->ibController.writeFluidVelocities();
  // logInfo("Map data");
  _im->ibController.mapData(_im->memory.timeStepSize());
  _im->ibController.readFluidForces();

  // logInfo("Read data");
  VectorDsType total_force = VectorDsType::Zero();

  for (auto it = _im->ibController.begin();
       it != _im->ibController.end();
       ++it) {
    auto force = _im->ibController.getFluidForce(it);

    _im->memory.fgh()[it->first] += force;

    auto body_force = _im->ibController.computeBodyForceAt(
      it->first,
      force,
      &_im->memory);

    total_force += body_force;
  }

  Private::mpiAllReduce<ScalarType>(MPI_IN_PLACE,
                                    total_force.data(),
                                    Dimensions,
                                    MPI_SUM,
                                    PETSC_COMM_WORLD);

  _im->reporter->addAt(3, total_force);

  // logInfo("Vpe");
  _im->peSolver.executeVpe();
  _im->ghostHandlers.executeFghMpiExchange();

  // logInfo("Ppe");
  _im->peSolver.executePpe();
  _im->ghostHandlers.executePressureMpiExchange();

  for (auto& accessor : _im->memory.grid()->innerGrid) {
    VectorDsType grad_pressure
      = PressureProcessing<SolverId>::grad(accessor);

    accessor.velocity()
      = accessor.fgh()
        - _im->memory.timeStepSize() * grad_pressure;

    compute_max_velocity(accessor, _im->memory.maxVelocity());
  }

  _im->ghostHandlers.executeVelocityInitialization();
  _im->ghostHandlers.executeVelocityMpiExchange();

  _im->ibController.computeBodyForce(&_im->memory, _im->reporter);

  _im->memory.time() += _im->memory.timeStepSize();
  ++_im->memory.iterationNumber();

  _im->reporter->addAt(0, _im->memory.iterationNumber());
  _im->reporter->addAt(1, _im->memory.time());
  _im->reporter->addAt(2, _im->memory.timeStepSize());
}

Fluid_InstantiateExternTemplates(FsfdSolver);

}
}
