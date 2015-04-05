#include "facade.hpp"

#include "Configuration.hpp"
#include "ParticularSimulationController.hpp"
#include "SolverBuilder.hpp"
#include "SolverBuilderTraits.hpp"
#include "VtkOutput/Writer.hpp"
#include "XdmfHdf5Output/Writer.hpp"

#include <Uni/ExecutionControl/exception>

//

namespace FsiSimulation {
namespace FluidSimulation {
namespace Private {
template <typename TSolverBuilderTraits, typename TSimulationController>
inline std::unique_ptr<SimulationController>
_create_simulation_controller(FluidSimulation::Configuration* configuration) {
  using namespace FluidSimulation;

  using SimulationControllerType = TSimulationController;
  using SolverBuilderTraitsType  = TSolverBuilderTraits;
  using SolverTraitsType
          = typename SimulationControllerType::SolverTraitsType;

  static_assert((SolverTraitsType::Dimensions > 1)
                && (SolverTraitsType::Dimensions < 4),
                "Only 2D and 3D simulations supported");

  std::unique_ptr<SimulationControllerType> controller(
    new SimulationControllerType());

  SolverBuilder<SolverBuilderTraitsType> builder(configuration,
                                                 controller->solver());

  controller->iterationLimit() = configuration->iterationLimit;
  controller->timeLimit()      = configuration->timeLimit;
  controller->plotInterval()   = configuration->plotInterval;

  return std::unique_ptr<SimulationController>(controller.release());
}

template <int TDimensions,
          typename TScalar,
          int TSolverType,
          int TImmersedBoudnaryType,
          int TDebug>
inline std::unique_ptr<SimulationController>
_create_simulation_controller_nd_scalar_solver_type_immersed_boudnary_type_debug(
  FluidSimulation::Configuration* configuration) {
  using SolverBuilderTraitsType
          = SolverBuilderTraits<TScalar,
                                TDimensions,
                                TSolverType,
                                TImmersedBoudnaryType,
                                TDebug>;

  using SolverTraitsType
          = typename SolverBuilderTraitsType::SolverTraitsType;

  using MemoryType
          = typename SolverTraitsType::MemoryType;

  // if (configuration->outputType == OutputEnum::Vtk) {
  //   using SimulationControllerType
  //           = ParticularSimulationController
  //             < SolverTraitsType, VtkOutput::Writer < MemoryType >>;

  //   return _create_simulation_controller<SolverBuilderTraitsType,
  //                                        SimulationControllerType>
  //            (configuration);
  // } else if (configuration->outputType == OutputEnum::Xdmf) {
    using SimulationControllerType
            = ParticularSimulationController
              < SolverTraitsType, XdmfHdf5Output::Writer < MemoryType >>;

    return _create_simulation_controller<SolverBuilderTraitsType,
                                         SimulationControllerType>
             (configuration);
  // }
  // throwException(
  //   "Failed to crate simulation controller for the provided output type");
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
     1>(configuration);
}

template <int TDimensions, typename TScalar, int TSolverType>
inline std::unique_ptr<SimulationController>
_create_simulation_controller_nd_scalar_solver_type(
  FluidSimulation::Configuration* configuration) {
  if (configuration->doImmersedBoundary) {
    return
      _create_simulation_controller_nd_scalar_solver_type_immersed_boudnary_type
      <TDimensions, TScalar, TSolverType, 1>(configuration);
  }

  return
    _create_simulation_controller_nd_scalar_solver_type_immersed_boudnary_type
    <TDimensions, TScalar, TSolverType, 1>(configuration);
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
           <TDimensions, TScalar, 0>(configuration);
  }
  throwException(
    "Failed to crate simulation controller for the provided solver type");
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
}
}

std::unique_ptr<SimulationController>
create_simulation_controller(FluidSimulation::Configuration* configuration) {
  using namespace FluidSimulation;

  if (configuration->dimensions == 2) {
    return Private::_create_simulation_controller_nd<2>(configuration);
  } else if (configuration->dimensions == 3) {
    return Private::_create_simulation_controller_nd<2>(configuration);
  }
  throwException(
    "Failed to create simulation controller for the provided dimension '{1}'"
    ", only 2 and 3 dimension numbers are allowed",
    configuration->dimensions);
}
}
}
