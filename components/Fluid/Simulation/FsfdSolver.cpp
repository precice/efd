#include "FsfdSolver.hpp"

#include "Configuration.hpp"
#include "GhostLayer/IfsfdHandlers.hpp"
#include "GhostLayer/SfsfdHandlers.hpp"
#include "Grid.hpp"
#include "GridGeometry.hpp"
#include "IfsfdCellAccessor.hpp"
#include "IfsfdMemory.hpp"
#include "ImmersedBoundary/BodyForce/functions.hpp"
#include "ImmersedBoundary/Controller.hpp"
#include "ImmersedBoundary/EmptyController.hpp"
#include "ImmersedBoundary/PreciceBasedController.hpp"
#include "ImmersedBoundary/RbfBasedController.hpp"
#include "ImmersedBoundary/functions.hpp"
#include "ImmersedBoundary/processcoupling.hpp"
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

#include <functional>
#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TVector>
inline typename TVector::Scalar
compute_cfl(TVector const& minCellWidth,
            TVector const& maxVelocity) {
  using Scalar = typename TVector::Scalar;

  Scalar localMin    = std::numeric_limits<Scalar>::max();
  bool   is_computed = false;

  for (unsigned d = 0; d < TVector::RowsAtCompileTime; ++d) {
    if (std::abs(maxVelocity(d)) > std::numeric_limits<Scalar>::epsilon()) {
      localMin    = std::min(localMin, 1.0 / maxVelocity(d));
      is_computed = true;
    }
  }

  if (!is_computed) {
    return minCellWidth.minCoeff();
  }

  return localMin;
}

template <typename TScalar, typename TVector>
inline TScalar
compute_explicit_time_step_size(TScalar const& re,
                                TScalar const& tau,
                                TVector const& minCellWidth,
                                TVector const& maxVelocity) {
  TScalar factor = minCellWidth.cwiseProduct(minCellWidth).cwiseInverse().sum();

  TScalar localMin = (TScalar)(re / (2.0 * factor));

  localMin = std::min(localMin, compute_cfl(minCellWidth, maxVelocity));

  TScalar globalMin = std::numeric_limits<TScalar>::max();
  Private::mpi_all_reduce(&localMin,
                          &globalMin,
                          1,
                          MPI_MIN,
                          PETSC_COMM_WORLD);

  factor  = globalMin;
  factor *= tau;

  return factor;
}

template <typename TVector>
inline typename TVector::Scalar
compute_implicit_time_step_size(typename TVector::Scalar const& tau,
                                TVector const&                  minCellWidth,
                                TVector const&                  maxVelocity) {
  using Scalar = typename TVector::Scalar;

  Scalar localMin = compute_cfl(minCellWidth, maxVelocity);

  Scalar globalMin = std::numeric_limits<Scalar>::max();
  Private::mpi_all_reduce(&localMin,
                          &globalMin,
                          1,
                          MPI_MIN,
                          PETSC_COMM_WORLD);

  globalMin *= tau;

  return globalMin;
}

template <typename TSolverTraits>
class FsfdSolverImplementation {
public:
  using Interface = FsfdSolver<TSolverTraits>;

  FsfdSolverImplementation(Interface* in)
    : _in(in),
    preciceInterface(nullptr),
    maxLayerSize(0),
    doCoupling(false) {}

  Uni_Firewall_INTERFACE_LINK(FsfdSolver<TSolverTraits> );

  std::function<void()> iterateFunction;
  std::function<void()> locateStructureFunction;
  typename Interface::MemoryType                     memory;
  typename Interface::PeSolverType                   peSolver;
  typename Interface::GhostHandlersType              ghostHandlers;
  precice::SolverInterface* preciceInterface;

  std::unique_ptr<typename Interface::IbControllerType> ibController;
  unsigned                                              maxLayerSize;
  bool                                                  doCoupling;

  Reporter* reporter;
};

template <typename T>
FsfdSolver<T>::
FsfdSolver(Configuration const* configuration) : _im(new Implementation(this)) {
  _im->memory.parameters()->re()
    = configuration->get<long double>("/Equations/Ins/ReynoldsNumber");

  if (!configuration->isOfType<std::string>(
        "/Equations/Ins/DiffusionMultiplier")) {
    _im->memory.parameters()->diffusionMultiplier()
      = configuration->get<long double>(
      "/Equations/Ins/DiffusionMultiplier");
  } else {
    _im->memory.parameters()->diffusionMultiplier()
      = 1.0 / _im->memory.parameters()->re();
  }

  if (!configuration->isOfType<std::string>(
        "/Equations/Ins/PressureGradientMultiplier")) {
    _im->memory.parameters()->gradPressureMultiplier()
      = configuration->get<long double>(
      "/Equations/Ins/PressureGradientMultiplier");
  } else {
    _im->memory.parameters()->gradPressureMultiplier()
      = 1.0;
  }

  if (configuration->is("/Ib/Features/FullVelocityPrediction")) {
    _im->iterateFunction
      = std::bind(&FsfdSolver::iterateWithFullIbVelocityPrediction, this);
  } else {
    _im->iterateFunction
      = std::bind(&FsfdSolver::iterateWithFastIbVelocityPrediction, this);
  }

  if (configuration->is("/Ib/Features/DevelopingStructure")) {
    _im->locateStructureFunction
      = std::bind(&FsfdSolver::locateStructure, this);
  } else {
    _im->locateStructureFunction
      = [this] () {
          this->locateStructure();
          _im->locateStructureFunction = [] () {};
        };
  }

  if (configuration->is("/Ib/Features/Coupling")) {
    _im->doCoupling = true;
  }

  if (configuration->is("/Ib/Schemes/DirectForcing/PreciceBased")) {
    using PreciceBasedControllerType
            = ImmersedBoundary::PreciceBasedController<SolverTraitsType>;
    using UniquePreciceBasedControllerType
            = std::unique_ptr<PreciceBasedControllerType>;

    UniquePreciceBasedControllerType temp(
      new PreciceBasedControllerType(configuration, &_im->memory));

    _im->maxLayerSize = std::max(_im->maxLayerSize, temp->getMaxLayerSize());

    _im->ibController.reset(temp.release());
  } else if (configuration->is("/Ib/Schemes/DirectForcing/RbfBased")) {
    _im->ibController.reset(
      new ImmersedBoundary::RbfBasedController<SolverTraitsType>(
        configuration,
        &_im->memory));
  } else {
    _im->ibController.reset(
      new ImmersedBoundary::EmptyController<VectorDsType> );
  }
}

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
initialize(precice::SolverInterface* precice_interface,
           Reporter*                 reporter) {
  _im->preciceInterface = precice_interface;

  _im->maxLayerSize = std::max(_im->maxLayerSize, 1u);

  if (!precice_interface) {
    _im->ibController.reset(
      new ImmersedBoundary::EmptyController<VectorDsType> );
  }

  _im->ibController->initialize(_im->preciceInterface);

  _im->reporter = reporter;

  _im->peSolver.initialize(&_im->memory, &_im->ghostHandlers);

  _im->memory.maxVelocity()     = VectorDsType::Zero();
  _im->memory.timeStepSize()    = 0.0;
  _im->memory.time()            = 0.0;
  _im->memory.iterationNumber() = 0;

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

  _im->reporter->addAt("IterationNumber", _im->memory.iterationNumber());
  _im->reporter->addAt("Time", _im->memory.time());
  _im->reporter->addAt("TimeStepSize", _im->memory.timeStepSize());
  _im->reporter->addAt("IbForceSum", VectorDsType::Zero().eval());
  _im->reporter->addAt("Force1", VectorDsType::Zero().eval());
  _im->reporter->addAt("Force2", VectorDsType::Zero().eval());
  _im->reporter->addAt("Force3", VectorDsType::Zero().eval());
  _im->reporter->addAt("MinVelocity", VectorDsType::Zero().eval());
  _im->reporter->addAt("MaxVelocity", VectorDsType::Zero().eval());
  _im->reporter->addAt("MinPressure", 0);
  _im->reporter->addAt("MaxPressure", 0);
  _im->reporter->recordIteration();
}

template <typename T>
typename FsfdSolver<T>::ScalarType
FsfdSolver<T>::
computeTimeStepSize() {
  if (SolverId == 0) {
    _im->memory.timeStepSize()
      = compute_explicit_time_step_size(_im->memory.parameters()->re(),
                                        _im->memory.parameters()->tau(),
                                        _im->memory.gridGeometry()->minCellWidth(),
                                        _im->memory.maxVelocity());
  } else {
    _im->memory.timeStepSize()
      = compute_implicit_time_step_size(_im->memory.parameters()->tau(),
                                        _im->memory.gridGeometry()->minCellWidth(),
                                        _im->memory.maxVelocity());
  }

  if (_im->memory.parallelDistribution()->rank() == 0) {
    logInfo("dt = {1}", _im->memory.timeStepSize());
    logInfo("maxv = {1}",
            _im->memory.maxVelocity().cwiseProduct(
              _im->memory.gridGeometry()->minCellWidth()).transpose());
  }
  _im->memory.maxVelocity()
    = VectorDsType::Constant(std::numeric_limits<ScalarType>::min());

  return _im->memory.timeStepSize();
}

template <typename T>
void
FsfdSolver<T>::
iterate() {
  // logInfo("Locate structure ...");
  _im->locateStructureFunction();
  // logInfo("Locate structure is finished");
  _im->iterateFunction();
}

template <typename T>
void
FsfdSolver<T>::
iterateWithFastIbVelocityPrediction() {
  for (auto const& accessor : _im->memory.grid()->innerGrid) {
    auto parts = updateFgh(accessor);

    auto ib_fluid_cell
      = _im->ibController->getVelocityIterable().find(accessor.globalIndex());

    if (ib_fluid_cell != _im->ibController->getVelocityIterable().end()) {
      auto temp
        = accessor.velocity()
          + _im->memory.timeStepSize() * (-1.5 * parts[0]
                                          + 0.5 * accessor.convection()
                                          + parts[1]
                                          - parts[2]
                                          + _im->memory.parameters()->g());

      ib_fluid_cell->data() = temp;
      // logInfo("Add velocity {1}", ib_fluid_cell->data().transpose());
    }

    accessor.convection() = parts[0];
  }

  advanceFsi();

  addIbForces();

  solvePoissonEquations();

  updateVelocities();

  finalizeIteration();
}

template <typename T>
void
FsfdSolver<T>::
iterateWithFullIbVelocityPrediction() {
  logInfo("Start Deep Ib velocity prediction ...");

  for (auto const& accessor : _im->memory.grid()->innerGrid) {
    updateFgh(accessor);
  }

  _im->peSolver.executeVpe();
  _im->ghostHandlers.executeFghMpiExchange();

  // std::unique_ptr<ScalarType[]> pressure_backup(
  // new ScalarType[_im->memory.grid()->size().prod()]);

  // std::memcpy(pressure_backup.get(),
  // _im->memory.pressure(),
  // _im->memory.grid()->size().prod() * sizeof (ScalarType));

  // solvePoissonEquations();

  for (auto& ib_fluid_cell : _im->ibController->getVelocityIterable()) {
    auto accessor = *_im->memory.grid()->begin();
    accessor.initialize(ib_fluid_cell.globalIndex());

    ib_fluid_cell.data() = accessor.fgh();

    // VectorDsType grad_pressure
    // = PressureProcessing<SolverId>::grad(accessor);

    // ib_fluid_cell.data()
    // = accessor.fgh()
    // - _im->memory.timeStepSize()
    // * _im->memory.parameters()->gradPressureMultiplier() * grad_pressure;
  }

  // std::memcpy(_im->memory.pressure(),
  // pressure_backup.get(),
  // _im->memory.grid()->size().prod() * sizeof (ScalarType));

  _im->memory.maxVelocity()
    = VectorDsType::Constant(std::numeric_limits<ScalarType>::min());

  logInfo("Finish the prediction phase");

  logInfo("Start solving phase ...");

  for (auto const& accessor : _im->memory.grid()->innerGrid) {
    auto parts = updateFgh(accessor);

    accessor.convection() = parts[0];
  }

  advanceFsi();

  addIbForces();

  solvePoissonEquations();

  updateVelocities();

  finalizeIteration();

  logInfo("Finish deep Ib velocity prediction");
}

template <typename T>
std::array<typename FsfdSolver<T>::VectorDsType, 3>
FsfdSolver<T>::
updateFgh(CellAccessorType const& accessor) {
  std::array<VectorDsType, 3> parts;

  _im->memory.setForceAt(accessor.globalIndex(), VectorDsType::Zero());

  parts[0] = ConvectionProcessing<Dimensions>::compute(
    accessor,
    _im->memory.parameters()->gamma());

  parts[1] = DiffusionProcessing<Dimensions>::compute(accessor);

  accessor.setDiffusion(parts[1]);

  parts[1] = _im->memory.parameters()->diffusionMultiplier() * parts[1];

  for (int d = 0; d < Dimensions; ++d) {
    parts[2](d)
      = (accessor.pressure(d, +1) - accessor.pressure())
        / (0.5 * (accessor.width(d, +1, d) + accessor.width(d)));
  }

  parts[2] = _im->memory.parameters()->gradPressureMultiplier() * parts[2];

  if (SolverId == 0) {
    accessor.fgh()
      = accessor.velocity()
        + _im->memory.timeStepSize() * (parts[1]
                                        - parts[0]
                                        + _im->memory.parameters()->g()
                                        );
  } else {
    accessor.fgh()
      = accessor.velocity()
        + _im->memory.timeStepSize() * (-1.5 * parts[0]
                                        + 0.5 * accessor.convection()
                                        + 0.5 * parts[1]
                                        - parts[2]
                                        + _im->memory.parameters()->g());
  }

  return parts;
}

template <typename T>
void
FsfdSolver<T>::
advanceFsi() {
  namespace pc = precice::constants;
  std::string writeCheckpoint(pc::actionWriteIterationCheckpoint());
  std::string readCheckpoint(pc::actionReadIterationCheckpoint());

  if (_im->preciceInterface != nullptr) {
    if (_im->preciceInterface->isActionRequired(writeCheckpoint)) {
      _im->preciceInterface->fulfilledAction(writeCheckpoint);
    }
  }

  sendCouplingData();

  _im->ibController->processVelocities();

  // logInfo("Invoking PreCICE's advance");
  if (_im->preciceInterface != nullptr) {
    _im->preciceInterface->advance(_im->memory.timeStepSize());
    // logInfo("Finished PreCICE's advance invokation");

    if (_im->preciceInterface->isActionRequired(readCheckpoint)) {
      _im->preciceInterface->fulfilledAction(readCheckpoint);
    }
  }
  // logInfo("Finished PreCICE's advance invokation");
  _im->ibController->processForces();
}

template <typename T>
void
FsfdSolver<T>::
addIbForces() {
  // logInfo("Read data");
  VectorDsType total_force = VectorDsType::Zero();

  for (auto& ib_fluid_cell : _im->ibController->getForceIterable()) {
    _im->memory.fgh()[ib_fluid_cell.globalIndex()]
      += _im->memory.timeStepSize() * ib_fluid_cell.data();

    // logInfo("Add force {1}", ib_fluid_cell.data().transpose());

    auto accessor = *_im->memory.grid()->begin();
    accessor.initialize(ib_fluid_cell.globalIndex());

    VectorDsType body_force = accessor.width().prod()
                              * (-ib_fluid_cell.data());

    _im->memory.setForceAt(ib_fluid_cell.globalIndex(), body_force);

    total_force += body_force;
  }
  Private::mpi_all_reduce(MPI_IN_PLACE,
                          total_force.data(),
                          Dimensions,
                          MPI_SUM,
                          PETSC_COMM_WORLD);

  _im->reporter->addAt("IbForceSum", total_force);
}

template <typename T>
void
FsfdSolver<T>::
solvePoissonEquations() {
  // logInfo("Vpe");
  _im->peSolver.executeVpe();
  _im->ghostHandlers.executeFghMpiExchange();

  // logInfo("Ppe");
  _im->peSolver.executePpe();
  _im->ghostHandlers.executePressureMpiExchange();
}

template <typename T>
void
FsfdSolver<T>::
updateVelocities() {
  VectorDsType min_velocity
    = VectorDsType::Constant(std::numeric_limits<ScalarType>::max());
  VectorDsType max_velocity
    = VectorDsType::Constant(std::numeric_limits<ScalarType>::min());

  ScalarType min_pressure = std::numeric_limits<ScalarType>::max();
  ScalarType max_pressure = std::numeric_limits<ScalarType>::min();

  for (auto& accessor : _im->memory.grid()->innerGrid) {
    VectorDsType grad_pressure
      = PressureProcessing<SolverId>::grad(accessor);

    accessor.velocity()
      = accessor.fgh()
        - _im->memory.timeStepSize()
        * _im->memory.parameters()->gradPressureMultiplier() * grad_pressure;

    compute_max_velocity(accessor, _im->memory.maxVelocity());

    if ((accessor.positionInRespectToGeometry().array() > 0).all()) {
      min_velocity = accessor.velocity().cwiseMin(min_velocity);
      max_velocity = accessor.velocity().cwiseMax(max_velocity);
      min_pressure = std::min(accessor.pressure(), min_pressure);
      max_pressure = std::max(accessor.pressure(), max_pressure);
    }
  }

  Private::mpi_all_reduce(MPI_IN_PLACE,
                          min_velocity.data(),
                          Dimensions,
                          MPI_MIN,
                          _im->memory.parallelDistribution()->mpiCommunicator);
  Private::mpi_all_reduce(MPI_IN_PLACE,
                          max_velocity.data(),
                          Dimensions,
                          MPI_MAX,
                          _im->memory.parallelDistribution()->mpiCommunicator);
  Private::mpi_all_reduce(MPI_IN_PLACE,
                          &min_pressure,
                          1,
                          MPI_MIN,
                          _im->memory.parallelDistribution()->mpiCommunicator);
  Private::mpi_all_reduce(MPI_IN_PLACE,
                          &max_pressure,
                          1,
                          MPI_MAX,
                          _im->memory.parallelDistribution()->mpiCommunicator);

  _im->ghostHandlers.executeVelocityInitialization();
  _im->ghostHandlers.executeVelocityMpiExchange();

  _im->reporter->addAt("MinVelocity", min_velocity);
  _im->reporter->addAt("MaxVelocity", max_velocity);
  _im->reporter->addAt("MinPressure", min_pressure);
  _im->reporter->addAt("MaxPressure", max_pressure);
}

template <typename T>
void
FsfdSolver<T>::
finalizeIteration() {
  computeBodyForce();

  _im->memory.time() += _im->memory.timeStepSize();
  ++_im->memory.iterationNumber();

  _im->reporter->addAt("IterationNumber", _im->memory.iterationNumber());
  _im->reporter->addAt("Time", _im->memory.time());
  _im->reporter->addAt("TimeStepSize", _im->memory.timeStepSize());
}

template <typename T>
void
FsfdSolver<T>::
locateStructure() {
  namespace ib = ImmersedBoundary;

  if (_im->preciceInterface == nullptr) {
    return;
  }

  if (!_im->preciceInterface->hasMesh("BodyMesh")) {
    throwException("Precice configuration does not have 'BodyMesh'");
  }

  std::set<int> mesh_set({ _im->preciceInterface->getMeshID("BodyMesh") });

  // logInfo("Start computing distance of cells from structure surface");

  for (auto const& accessor : * _im->memory.grid()) {
    // logInfo("Compute distance for one cell is being started ...");
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

  this->_im->ghostHandlers.executeLocationsMpiExchange();

  // logInfo("End computing distance of cells from structure surface");

  for (auto const& accessor : * _im->memory.grid()) {
    ib::set_cell_neighbors_along_geometry_interface(
      accessor,
      _im->preciceInterface,
      mesh_set,
      _im->maxLayerSize);
  }
  _im->ibController->precompute();

  // Coupling mesh construction

  if (!_im->doCoupling) {
    return;
  }
  namespace ib = ImmersedBoundary;

  ib::create_coupling_mesh(_im->memory.grid()->innerGrid,
                           _im->preciceInterface);
}

template <typename T>
void
FsfdSolver<T>::
sendCouplingData() {
  if (_im->preciceInterface == nullptr) {
    return;
  }

  if (!_im->doCoupling) {
    return;
  }
  namespace ib = ImmersedBoundary;

  ib::send_coupling_stresses(_im->memory.grid()->innerGrid,
                             _im->preciceInterface,
                             _im->memory.parameters()->diffusionMultiplier(),
                             _im->memory.parameters()->gradPressureMultiplier());
}

template <typename T>
void
FsfdSolver<T>::
computeBodyForce() {
  VectorDsType total_force       = VectorDsType::Zero();
  VectorDsType total_force_v2    = VectorDsType::Zero();
  VectorDsType total_force_turek = VectorDsType::Zero();

  for (auto const& accessor : _im->memory.grid()->innerGrid) {
    namespace fc = ImmersedBoundary::BodyForce;
    VectorDsType force       = VectorDsType::Zero();
    VectorDsType force_v2    = VectorDsType::Zero();
    VectorDsType force_turek = VectorDsType::Zero();
    fc::compute_cell_force(
      accessor,
      _im->memory.parameters()->diffusionMultiplier(),
      _im->memory.parameters()->gradPressureMultiplier(),
      force);
    fc::compute_cell_force_v2(
      accessor,
      _im->memory.parameters()->diffusionMultiplier(),
      _im->memory.parameters()->gradPressureMultiplier(),
      force_v2);
    fc::compute_cell_force_turek(
      accessor,
      _im->memory.parameters()->diffusionMultiplier(),
      _im->memory.parameters()->gradPressureMultiplier(),
      force_turek);
    accessor.setBodyForce(force);
    total_force       += force;
    total_force_v2    += force_v2;
    total_force_turek += force_turek;
  }

  FluidSimulation::Private::mpi_all_reduce(MPI_IN_PLACE,
                                           total_force.data(),
                                           Dimensions,
                                           MPI_SUM,
                                           PETSC_COMM_WORLD);

  FluidSimulation::Private::mpi_all_reduce(MPI_IN_PLACE,
                                           total_force_v2.data(),
                                           Dimensions,
                                           MPI_SUM,
                                           PETSC_COMM_WORLD);

  FluidSimulation::Private::mpi_all_reduce(MPI_IN_PLACE,
                                           total_force_turek.data(),
                                           Dimensions,
                                           MPI_SUM,
                                           PETSC_COMM_WORLD);

  _im->reporter->addAt("Force1", total_force);
  _im->reporter->addAt("Force2", total_force_v2);
  _im->reporter->addAt("Force3", total_force_turek);
}

Fluid_InstantiateExternTemplates(FsfdSolver);
}
}
