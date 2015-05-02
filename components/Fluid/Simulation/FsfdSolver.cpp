#include "FsfdSolver.hpp"

#include "GhostLayer/IfsfdHandlers.hpp"
#include "GhostLayer/SfsfdHandlers.hpp"
#include "Grid.hpp"
#include "GridGeometry.hpp"
#include "IfsfdCellAccessor.hpp"
#include "IfsfdMemory.hpp"
#include "ImmersedBoundary/BodyForce/functions.hpp"
#include "ImmersedBoundary/Controller.hpp"
#include "ImmersedBoundary/PreciceBasedController.hpp"
#include "ImmersedBoundary/RbfBasedController.hpp"
#include "ImmersedBoundary/functions.hpp"
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
  precice::SolverInterface* preciceInterface;

  std::unique_ptr<typename Interface::IbControllerType> ibController;
  unsigned                                              maxLayerSize;

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
typename FsfdSolver<T>::IbControllerType const*
FsfdSolver<T>::
immersedBoundaryController() const {
  return _im->ibController.get();
}

template <typename T>
typename FsfdSolver<T>::IbControllerType *
FsfdSolver<T>::
immersedBoundaryController() {
  return _im->ibController.get();
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
initialize(precice::SolverInterface* preciceInterface,
           Reporter*                 reporter) {
  _im->preciceInterface = preciceInterface;
  _im->preciceInterface->initialize();
  _im->preciceInterface->initializeData();

  _im->maxLayerSize = 1;
  _im->ibController.reset(
    new ImmersedBoundary::RbfBasedController<SolverTraitsType>(
      _im->preciceInterface,
      &_im->memory));
  _im->ibController->initialize();

  _im->reporter = reporter;

  _im->peSolver.initialize(&_im->memory, &_im->ghostHandlers);

  _im->memory.maxVelocity()     = VectorDsType::Zero();
  _im->memory.timeStepSize()    = 0.0;
  _im->memory.time()            = 0.0;
  _im->memory.iterationNumber() = 0;

  locateInterfaceCells();
  _im->ibController->precompute();

  for (auto& accessor : _im->memory.grid()->innerGrid) {
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
    auto convection = ConvectionProcessing<Dimensions>::compute(
      accessor, _im->memory.parameters()->gamma());
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
                       1);
  // _im->ibController.outerLayerSize());
  _im->reporter->setAt("InnerLayerSize",
                       1);
  // _im->ibController.innerLayerSize());
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

    auto ib_fluid_cell
      = _im->ibController->getVelocityIterable().find(accessor.globalIndex());

    if (ib_fluid_cell != _im->ibController->getVelocityIterable().end()) {
      _im->memory.setForceAt(accessor.globalIndex(), VectorDsType::Zero());
      auto temp
        = accessor.velocity()
          + _im->memory.timeStepSize() * (-1.5 * convection
                                          + 0.5 * previousConvection
                                          + diffusion
                                          - grad_pressure
                                          + _im->memory.parameters()->g());

      ib_fluid_cell->data() = temp;
    }
  }

  _im->ibController->processVelocities();
  // logInfo("Map data");
  _im->preciceInterface->advance(_im->memory.timeStepSize());
  _im->ibController->processForces();

  // logInfo("Read data");
  VectorDsType total_force = VectorDsType::Zero();

  for (auto& ib_fluid_cell : _im->ibController->getForceIterable()) {
    _im->memory.fgh()[ib_fluid_cell.globalIndex()]
      += _im->memory.timeStepSize() * ib_fluid_cell.data();

    auto accessor = *_im->memory.grid()->begin();
    accessor.initialize(ib_fluid_cell.globalIndex());

    VectorDsType body_force = accessor.width().prod()
                              * ib_fluid_cell.data();

    _im->memory.setForceAt(ib_fluid_cell.globalIndex(), body_force);

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

  computeBodyForce();

  _im->memory.time() += _im->memory.timeStepSize();
  ++_im->memory.iterationNumber();

  _im->reporter->addAt(0, _im->memory.iterationNumber());
  _im->reporter->addAt(1, _im->memory.time());
  _im->reporter->addAt(2, _im->memory.timeStepSize());
}

template <typename T>
void
FsfdSolver<T>::
locateInterfaceCells() {
  namespace ib = ImmersedBoundary;

  if (!_im->preciceInterface->hasMesh("BodyMesh")) {
    throwException("Precice configuration does not have 'BodyMesh'");
  }

  std::set<int> mesh_set({ _im->preciceInterface->getMeshID("BodyMesh") });

  for (auto const& accessor : * _im->memory.grid()) {
    for (unsigned d = 0; d < Dimensions; ++d) {
      accessor.positionInRespectToGeometry(d)
        = ib::convert_precice_position(
        _im->preciceInterface->inquirePosition(
          accessor.velocityPosition(d).template cast<double>().data(),
          mesh_set));
    }
    accessor.positionInRespectToGeometry(Dimensions)
      = ib::convert_precice_position(
      _im->preciceInterface->inquirePosition(
        accessor.pressurePosition().template cast<double>().data(),
        mesh_set));
  }

  for (auto const& accessor : * _im->memory.grid()) {
    ib::set_cell_neighbors_along_geometry_interface(
      accessor,
      _im->preciceInterface,
      mesh_set,
      _im->maxLayerSize);
  }
}

template <typename T>
void
FsfdSolver<T>::
computeBodyForce() {
  VectorDsType total_force       = VectorDsType::Zero();
  VectorDsType total_force_turek = VectorDsType::Zero();

  for (auto const& accessor : _im->memory.grid()->innerGrid) {
    namespace fc = ImmersedBoundary::BodyForce;
    VectorDsType force       = VectorDsType::Zero();
    VectorDsType force_turek = VectorDsType::Zero();
    fc::compute_cell_force(accessor,
                           _im->memory.parameters()->re(),
                           force);
    fc::compute_cell_force_turek(accessor,
                                 _im->memory.parameters()->re(),
                                 force_turek);
    accessor.setBodyForce(force);
    total_force       += force;
    total_force_turek += force_turek;
  }

  FluidSimulation::Private::mpiAllReduce<ScalarType>(MPI_IN_PLACE,
                                                     total_force.data(),
                                                     Dimensions,
                                                     MPI_SUM,
                                                     PETSC_COMM_WORLD);

  FluidSimulation::Private::mpiAllReduce<ScalarType>(MPI_IN_PLACE,
                                                     total_force_turek.data(),
                                                     Dimensions,
                                                     MPI_SUM,
                                                     PETSC_COMM_WORLD);

  _im->reporter->addAt(4, total_force);
  _im->reporter->addAt(5, total_force_turek);
}

Fluid_InstantiateExternTemplates(FsfdSolver);
}
}
