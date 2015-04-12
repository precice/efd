#include "FsfdSolver.hpp"

#include "Private/mpigenerics.hpp"
#include "Reporter.hpp"
#include "computemaxvelocity.hpp"
#include "functions.hpp"

#include "Grid.hpp"
#include "GridGeometry.hpp"
#include "IfsfdCellAccessor.hpp"
#include "SfsfdCellAccessor.hpp"

#include <precice/SolverInterface.hpp>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename T>
void
FsfdSolver<T>::
initialize(precice::SolverInterface* preciceInteface,
           Reporter*                 reporter) {
  _ibController.initialize(preciceInteface);
  _reporter = reporter;

  _peSolver.initialize(&_memory, &_ghostHandlers);

  _memory.maxVelocity()     = VectorDsType::Zero();
  _memory.timeStepSize()    = 0.0;
  _memory.time()            = 0.0;
  _memory.iterationNumber() = 0;

  for (auto& accessor : _memory.grid()->innerGrid) {
    _ibController.computePositionInRespectToGeometry(accessor);
    accessor.velocity() = VectorDsType::Zero();
    accessor.setDiffusion(VectorDsType::Zero());
    accessor.setForce(VectorDsType::Zero());
    accessor.setBodyForce(VectorDsType::Zero());
    PressureProcessing<SolverId>::initialize(accessor, 0.0);
    compute_max_velocity(accessor, _memory.maxVelocity());
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

  _reporter->setAt("Re",
                   _memory.parameters()->re());
  _reporter->setAt("Gamma",
                   _memory.parameters()->gamma());
  _reporter->setAt("Tau",
                   _memory.parameters()->tau());
  _reporter->setAt("G",
                   _memory.parameters()->g());
  _reporter->setAt("OuterLayerSize",
                   _ibController.outerLayerSize());
  _reporter->setAt("InnerLayerSize",
                   _ibController.innerLayerSize());
  _reporter->setAt("ProcessorSize",
                   _memory.parallelDistribution()->processorSize);
  _reporter->setAt("GlobalCellSize",
                   _memory.parallelDistribution()->globalCellSize);
  _reporter->setAt("UniformLocalCellSize",
                   _memory.parallelDistribution()->uniformLocalCellSize);
  _reporter->setAt("LastLocalCellSize",
                   _memory.parallelDistribution()->lastLocalCellSize);
  _reporter->setAt("Width",
                   _memory.gridGeometry()->size());
  _reporter->setAt("CellWidth",
                   _memory.gridGeometry()->minCellWidth());
  _reporter->setAt("SolverId", SolverId);

  _reporter->addAt(0, "IterationNumber");
  _reporter->addAt(1, "Time");
  _reporter->addAt(2, "TimeStepSize");
  _reporter->addAt(3, "Force1");
  _reporter->addAt(4, "Force2");
  _reporter->addAt(5, "Force3");

  _reporter->recordIteration();

  _reporter->addAt(0, _memory.iterationNumber());
  _reporter->addAt(1, _memory.time());
  _reporter->addAt(2, _memory.timeStepSize());
  _reporter->addAt(3, VectorDsType::Zero().eval());
  _reporter->addAt(4, VectorDsType::Zero().eval());
  _reporter->addAt(5, VectorDsType::Zero().eval());

  _reporter->recordIteration();
}

template <typename T>
void
FsfdSolver<T>::
iterate() {
  _memory.timeStepSize()
    = compute_time_step_size(_memory.parameters()->re(),
                             _memory.parameters()->tau(),
                             _memory.gridGeometry()->minCellWidth(),
                             _memory.maxVelocity());

  logInfo("dt = {1}", _memory.timeStepSize());
  logInfo("maxv = {1}",
          _memory.maxVelocity().cwiseProduct(_memory.gridGeometry()->minCellWidth()).transpose());
  _memory.maxVelocity()
    = VectorDsType::Constant(std::numeric_limits<ScalarType>::min());

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

    if (SolverId == 0) {
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
      // logInfo("Fgh | {1}", accessor.index().transpose());
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

  // logInfo("Map data");
  _ibController.mapData(_memory.timeStepSize());

  // logInfo("Read data");
  VectorDsType total_force = VectorDsType::Zero();

  for (auto it = _ibController.begin(); it != _ibController.end(); ++it) {
    VectorDsType force;
    _ibController.readFluidForce(it, force);

    _memory.fgh()[it->first] += force;

    auto body_force = _ibController.computeBodyForceAt(
      it->first,
      force,
      &_memory);

    total_force += body_force;
  }

  Private::mpiAllReduce<ScalarType>(MPI_IN_PLACE,
                                    total_force.data(),
                                    Dimensions,
                                    MPI_SUM,
                                    PETSC_COMM_WORLD);

  _reporter->addAt(3, total_force);

  // logInfo("Vpe");
  _peSolver.executeVpe();
  _ghostHandlers.executeFghMpiExchange();

  // logInfo("Ppe");
  _peSolver.executePpe();
  _ghostHandlers.executePressureMpiExchange();

  for (auto& accessor : _memory.grid()->innerGrid) {
    VectorDsType grad_pressure
      = PressureProcessing<SolverId>::grad(accessor);

    accessor.velocity()
      = accessor.fgh()
        - _memory.timeStepSize() * grad_pressure;

    compute_max_velocity(accessor, _memory.maxVelocity());
  }

  _ghostHandlers.executeVelocityInitialization();
  _ghostHandlers.executeVelocityMpiExchange();

  _ibController.computeBodyForce(&_memory, _reporter);

  _memory.time() += _memory.timeStepSize();
  ++_memory.iterationNumber();

  _reporter->addAt(0, _memory.iterationNumber());
  _reporter->addAt(1, _memory.time());
  _reporter->addAt(2, _memory.timeStepSize());
}

template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 0, double, 2 >>;
template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 1, double, 2 >>;
template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 0, double, 2 >>;
template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 1, double, 2 >>;
template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 0, double, 3 >>;
template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 1, double, 3 >>;
template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 0, double, 3 >>;
template class FsfdSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 1, double, 3 >>;

template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 0, double, 2 >>;
template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 1, double, 2 >>;
template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 0, double, 2 >>;
template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 1, double, 2 >>;
template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 0, double, 3 >>;
template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 1, double, 3 >>;
template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 0, double, 3 >>;
template class FsfdSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 1, double, 3 >>;
}
}
