#include "facade.hpp"

#include "Configuration.hpp"
#include "ParticularSimulationController.hpp"
#include "SolverBuilder.hpp"
#include "SolverBuilderTraits.hpp"
#include "VtkOutput/createinctance.hpp"
#include "XdmfHdf5Output/createinctance.hpp"
#include "SfsfdMemory.hpp"
#include "IfsfdMemory.hpp"
#include "FsfdSolver.hpp"

#include <Uni/ExecutionControl/exception>

namespace FsiSimulation {
namespace FluidSimulation {
namespace Private {
template <int TDimensions,
          typename TScalar,
          int TSolverType,
          int TImmersedBoudnaryType,
          int TDebug>
inline std::unique_ptr<SimulationController>
_create_simulation_controller_nd_scalar_solver_type_immersed_boudnary_type_debug(
  FluidSimulation::Configuration* configuration) {
  using SolverBuilderTraitsType
          = SolverBuilderTraits<TSolverType,
                                TImmersedBoudnaryType,
                                TDebug,
                                TScalar,
                                TDimensions>;

  using SolverTraitsType
          = typename SolverBuilderTraitsType::SolverTraitsType;

  using SimulationControllerType
          = ParticularSimulationController<SolverTraitsType>;

  using namespace FluidSimulation;

  static_assert((SolverTraitsType::Dimensions > 1)
                && (SolverTraitsType::Dimensions < 4),
                "Only 2D and 3D simulations supported");

  std::unique_ptr<SimulationControllerType> controller(
    new SimulationControllerType());

  SolverBuilder<SolverBuilderTraitsType> builder(configuration,
                                                 controller->solver());

  std::unique_ptr<IterationResultWriter> iteration_result_writer;

  if (configuration->outputType == OutputEnum::Vtk) {
    iteration_result_writer
      = VtkOutput::create_instance(controller->solver()->memory());
  } else if (configuration->outputType == OutputEnum::Xdmf) {
    iteration_result_writer
      = XdmfHdf5Output::create_instance(controller->solver()->memory());
  }

  controller->setIterationResultWriter(iteration_result_writer.release());

  controller->iterationLimit() = configuration->iterationLimit;
  controller->timeLimit()      = configuration->timeLimit;
  controller->plotInterval()   = configuration->plotInterval;

  return std::unique_ptr<SimulationController>(controller.release());
}

template <int TDimensions,
          typename TScalar,
          int TSolverType,
          int TImmersedBoudnaryType>
inline std::unique_ptr<SimulationController>
_create_simulation_controller_nd_scalar_solver_type_immersed_boudnary_type(
  FluidSimulation::Configuration* configuration) {
  if (configuration->doDebug) {
    return
      _create_simulation_controller_nd_scalar_solver_type_immersed_boudnary_type_debug
      <TDimensions,
       TScalar,
       TSolverType,
       TImmersedBoudnaryType,
       1>(configuration);
  }

  return
    _create_simulation_controller_nd_scalar_solver_type_immersed_boudnary_type_debug
    <TDimensions,
     TScalar,
     TSolverType,
     TImmersedBoudnaryType,
     0>(configuration);
}

template <int TDimensions, typename TScalar, int TSolverType>
inline std::unique_ptr<SimulationController>
_create_simulation_controller_nd_scalar_solver_type(
  FluidSimulation::Configuration* configuration) {
  if (configuration->doImmersedBoundary) {
    return
      _create_simulation_controller_nd_scalar_solver_type_immersed_boudnary_type
      <TDimensions, TScalar, TSolverType, 2>(configuration);
  }

  return
    _create_simulation_controller_nd_scalar_solver_type_immersed_boudnary_type
    <TDimensions, TScalar, TSolverType, 0>(configuration);
}

template <int TDimensions, typename TScalar>
inline std::unique_ptr<SimulationController>
_create_simulation_controller_nd_scalar(
  FluidSimulation::Configuration* configuration) {
  if (configuration->solverType == SolverEnum::Sfsfd) {
    return _create_simulation_controller_nd_scalar_solver_type
           <TDimensions, TScalar, 0>(configuration);
  } else if (configuration->solverType == SolverEnum::Ifsfd) {
    return _create_simulation_controller_nd_scalar_solver_type
           <TDimensions, TScalar, 1>(configuration);
  }
  throwException(
    "Failed to crate simulation controller for the provided solver type");

  return std::unique_ptr<SimulationController>();
}

template <int TDimensions>
inline std::unique_ptr<SimulationController>
_create_simulation_controller_nd(
  FluidSimulation::Configuration* configuration) {
  if (configuration->scalarType == ScalarEnum::Float) {
    return _create_simulation_controller_nd_scalar<TDimensions, double>(
      configuration);
  } else if (configuration->scalarType == ScalarEnum::Double) {
    return _create_simulation_controller_nd_scalar<TDimensions, double>(
      configuration);
  } else if (configuration->scalarType == ScalarEnum::LongDouble) {
    return _create_simulation_controller_nd_scalar<TDimensions, double>(
      configuration);
  }
  throwException(
    "Failed to create simulation controller for the provided configurations");

  return std::unique_ptr<SimulationController>();
}
}

std::unique_ptr<SimulationController>
create_simulation_controller(FluidSimulation::Configuration* configuration) {
  using namespace FluidSimulation;

  if (configuration->dimensions == 2) {
    return Private::_create_simulation_controller_nd<2>(configuration);
  } else if (configuration->dimensions == 3) {
    return Private::_create_simulation_controller_nd<3>(configuration);
  }
  throwException(
    "Failed to create simulation controller for the provided dimension '{1}'"
    ", only 2 and 3 dimension numbers are allowed",
    configuration->dimensions);

  return std::unique_ptr<SimulationController>();
}
}
}
