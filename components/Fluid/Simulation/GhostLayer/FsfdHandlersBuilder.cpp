#include "FsfdHandlersBuilder.hpp"

#include "InitializationActions.hpp"
#include "InitializationHandler.hpp"
#include "MpiExchangeHandler.hpp"
#include "PetscExchangeActions.hpp"
#include "PetscExchangeHandler.hpp"
#include "PressureStencilHanler.hpp"

#include "SfsfdHandlersBuilder.hpp"
#include "IfsfdHandlersBuilder.hpp"

#include "Simulation/Configuration.hpp"
#include "Simulation/FsfdSolver.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace Private {
template <typename T>
struct PpeTraits {
  typedef
    LsStencilGenerator::Handler
    <typename T::SolverTraitsType::GridType,
     T::Dimension,
     T::Direction>
    StencilGenerationHandler;

  using RhsGenerationAction
          = PetscExchange::PpeConstantRhsGenerationAction<
          typename T::SolverTraitsType::ScalarType>;
  typedef
    PetscExchange::Handler
    <typename T::SolverTraitsType::GridType,
     RhsGenerationAction,
     T::Dimension,
     T::Direction>
    RhsGenerationHandler;

  typedef
    typename T::PpeRhsAcquiererAction
    RhsAcquiererAction;
  typedef
    PetscExchange::Handler
    <typename T::SolverTraitsType::GridType,
     RhsAcquiererAction,
     T::Dimension,
     T::Direction>
    RhsAcquiererHandler;
};
}
template <typename T>
FsfdHandlersBuilder<T>::
FsfdHandlersBuilder(Configuration* configuration,
                    SolverType*    simulation) :
  _configuration(configuration),
  _solver(simulation) {}

template <typename T>
void
FsfdHandlersBuilder<T>::
setAsInput() {
  typedef
    Initialization::InputFghAction
    <SolverTraitsType,
     Dimension,
     Direction>
    InputFghInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType,
     InputFghInitializationAction,
     Dimensions,
     Dimension,
     Direction>
    InputFghInitialization;
  typedef
    Initialization::InputVelocityAction
    <SolverTraitsType, Dimension, Direction>
    InputVelocityInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType, InputVelocityInitializationAction, Dimensions,
     Dimension, Direction>
    InputVelocityInitialization;
  using Ppe = Private::PpeTraits<T>;

  _solver->ghostHandlers()->fghInitialization[Dimension][Direction]
    = InputFghInitialization::getHandler(
    &_solver->memory()->grid()
    ->indentedBoundaries[Dimension][Direction],
    _solver->memory()->parallelDistribution(),
    new InputFghInitializationAction(_configuration));

  _solver->ghostHandlers()->ppeStencilGeneratorStack[Dimension][Direction]
    = Ppe::StencilGenerationHandler::getNeumannMiddle(
    _solver->memory()->grid(),
    _solver->memory()->parallelDistribution());
  _solver->ghostHandlers()->ppeRhsGeneratorStack[Dimension][Direction]
    = Ppe::RhsGenerationHandler::getHandler(
    _solver->memory()->grid(),
    _solver->memory()->parallelDistribution(),
    new typename Ppe::RhsGenerationAction(0.0));

  if (Direction == 1) {
    _solver->ghostHandlers()->ppeRhsAcquiererStack[Dimension][Direction]
      = Ppe::RhsAcquiererHandler::getHandler(
      _solver->memory()->grid(),
      _solver->memory()->parallelDistribution(),
      new typename Ppe::RhsAcquiererAction());
  }

  _solver->ghostHandlers()->velocityInitialization[Dimension][Direction]
    = InputVelocityInitialization::getHandler(
    &_solver->memory()->grid()
    ->indentedBoundaries[Dimension][Direction],
    _solver->memory()->parallelDistribution(),
    new InputVelocityInitializationAction(
      _configuration,
      &_solver->memory()->maxVelocity()));
}

template <typename T>
void
FsfdHandlersBuilder<T>::
setAsParabolicInput() {
  typedef
    Initialization::ParabolicInputFghAction
    <SolverTraitsType,
     Dimension,
     Direction>
    ParabolicInputFghInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType,
     ParabolicInputFghInitializationAction,
     Dimensions,
     Dimension,
     Direction>
    ParabolicInputFghInitialization;
  typedef
    Initialization::ParabolicInputVelocityAction
    <SolverTraitsType, Dimension,
     Direction>
    ParabolicInputVelocityInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType, ParabolicInputVelocityInitializationAction,
     Dimensions,
     Dimension, Direction>
    ParabolicInputVelocityInitialization;
  using Ppe = Private::PpeTraits<T>;

  _solver->ghostHandlers()->fghInitialization[Dimension][Direction]
    = ParabolicInputFghInitialization::getHandler(
    &_solver->memory()->grid()
    ->indentedBoundaries[Dimension][Direction],
    _solver->memory()->parallelDistribution(),
    new ParabolicInputFghInitializationAction(_configuration));

  _solver->ghostHandlers()->ppeStencilGeneratorStack[Dimension][Direction]
    = Ppe::StencilGenerationHandler::getNeumannMiddle(
    _solver->memory()->grid(),
    _solver->memory()->parallelDistribution());
  _solver->ghostHandlers()->ppeRhsGeneratorStack[Dimension][Direction]
    = Ppe::RhsGenerationHandler::getHandler(
    _solver->memory()->grid(),
    _solver->memory()->parallelDistribution(),
    new typename Ppe::RhsGenerationAction(0.0));

  if (Direction == 1) {
    _solver->ghostHandlers()->ppeRhsAcquiererStack[Dimension][Direction]
      = Ppe::RhsAcquiererHandler::getHandler(
      _solver->memory()->grid(),
      _solver->memory()->parallelDistribution(),
      new typename Ppe::RhsAcquiererAction());
  }

  _solver->ghostHandlers()->velocityInitialization[Dimension][Direction]
    = ParabolicInputVelocityInitialization::getHandler(
    &_solver->memory()->grid()
    ->indentedBoundaries[Dimension][Direction],
    _solver->memory()->parallelDistribution(),
    new ParabolicInputVelocityInitializationAction(
      _configuration,
      &_solver->memory()->maxVelocity()));
}

template <typename T>
void
FsfdHandlersBuilder<T>::
setAsOutput() {
  typedef
    Initialization::OutputFghAction
    <SolverTraitsType,
     Dimension,
     Direction>
    OutputFghInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType,
     OutputFghInitializationAction,
     Dimensions,
     Dimension,
     Direction>
    OutputFghInitialization;
  typedef
    Initialization::OutputVelocityAction
    <SolverTraitsType, Dimension,
     Direction>
    OutputVelocityInitializationAction;
  typedef
    Initialization::Handler
    <typename GridType::BaseType, OutputVelocityInitializationAction,
     Dimensions,
     Dimension, Direction>
    OutputVelocityInitialization;
  using Ppe = Private::PpeTraits<T>;

  _solver->ghostHandlers()->fghInitialization[Dimension][Direction]
    = OutputFghInitialization::getHandler(
    &_solver->memory()->grid()
    ->indentedBoundaries[Dimension][Direction],
    _solver->memory()->parallelDistribution(),
    new OutputFghInitializationAction());

  _solver->ghostHandlers()->ppeStencilGeneratorStack[Dimension][Direction]
    = Ppe::StencilGenerationHandler::getDirichletMiddle(
    _solver->memory()->grid(),
    _solver->memory()->parallelDistribution());
  _solver->ghostHandlers()->ppeRhsGeneratorStack[Dimension][Direction]
    = Ppe::RhsGenerationHandler::getHandler(
    _solver->memory()->grid(),
    _solver->memory()->parallelDistribution(),
    new typename Ppe::RhsGenerationAction(0.0));

  if (Direction == 1) {
    _solver->ghostHandlers()->ppeRhsAcquiererStack[Dimension][Direction]
      = Ppe::RhsAcquiererHandler::getHandler(
      _solver->memory()->grid(),
      _solver->memory()->parallelDistribution(),
      new typename Ppe::RhsAcquiererAction());
  }

  _solver->ghostHandlers()->velocityInitialization[Dimension][Direction]
    = OutputVelocityInitialization::getHandler(
    &_solver->memory()->grid()
    ->indentedBoundaries[Dimension][Direction],
    _solver->memory()->parallelDistribution(),
    new OutputVelocityInitializationAction(
      &_solver->memory()->maxVelocity()));
}

template <typename T>
void
FsfdHandlersBuilder<T>::
setAsMpiExchange() {
  typedef
    MpiExchange::Handler
    <ScalarType,
     1,
     typename GridType::BaseType,
     FsfdHandlersBuilder::template getFgh<Dimension>,
     FsfdHandlersBuilder::template setFgh<Dimension>,
     Dimension,
     Direction>
    FghMpiExchangeHandler;
  typedef
    typename T::PressureMpiExchangeHandler
    PressureMpiExchangeHandler;
  typedef
    MpiExchange::Handler
    <ScalarType,
     Dimensions,
     typename GridType::BaseType,
     FsfdHandlersBuilder::getVelocity,
     FsfdHandlersBuilder::setVelocity,
     Dimension,
     Direction>
    VelocityMpiExchangeHandler;

  if (Direction == 0) {
    _solver->ghostHandlers()->mpiFghExchangeStack[Dimension][Direction]
      = FghMpiExchangeHandler::getReceiveHandler(
      &_solver->memory()->grid()->innerGrid,
      _solver->memory()->parallelDistribution());

    _solver->ghostHandlers()->mpiPressureExchangeStack[Dimension][Direction]
      = PressureMpiExchangeHandler::getSendHandler(
      &_solver->memory()->grid()->innerGrid,
      _solver->memory()->parallelDistribution());
  } else {
    _solver->ghostHandlers()->mpiFghExchangeStack[Dimension][Direction]
      = FghMpiExchangeHandler::getSendHandler(
      &_solver->memory()->grid()->innerGrid,
      _solver->memory()->parallelDistribution());

    _solver->ghostHandlers()->mpiPressureExchangeStack[Dimension][Direction]
      = PressureMpiExchangeHandler::getReceiveHandler(
      &_solver->memory()->grid()->innerGrid,
      _solver->memory()->parallelDistribution());
  }

  _solver->ghostHandlers()->mpiVelocityExchangeStack[Dimension][Direction]
    = VelocityMpiExchangeHandler::getExchangeHandler(
    &_solver->memory()->grid()->innerGrid,
    _solver->memory()->parallelDistribution());
}

template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 1 >>;

template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 1 >>;

template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 1 >>;

template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 1 >>;

template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 1 >>;

template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 1 >>;

template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 1 >>;

template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 1 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 0 >>;
template class FsfdHandlersBuilder
  < SfsfdHandlersBuilderTraits
  < SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 1 >>;

template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 1 >>;

template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 1 >>;

template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 1 >>;

template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 1 >>;

template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 1 >>;

template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 1 >>;

template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 1 >>;

template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 1 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 0 >>;
template class FsfdHandlersBuilder
  < IfsfdHandlersBuilderTraits
  < IfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 1 >>;
}
}
}
