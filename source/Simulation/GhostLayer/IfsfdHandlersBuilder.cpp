#include "IfsfdHandlersBuilder.hpp"

#include "PetscExchangeActions.hpp"
#include "PetscExchangeHandler.hpp"
#include "PressureStencilHanler.hpp"

#include "IfsfdHandlers.hpp"

#include "Simulation/Configuration.hpp"
#include "Simulation/FsfdSolver.hpp"
#include "Simulation/IfsfdMemory.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace Private {
template <typename SolverTraitsType,
          int TDimension,
          int TDirection>
struct VpeTraits {
  using GridType = typename SolverTraitsType::GridType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  typedef
    LsStencilGenerator::Handler
    <GridType, TDimension, TDirection>
    StencilGenerationHandler;

  using VxpeRhsAcquiererAction
          = PetscExchange::VpeRhsAcquiererAction<0, TDimension, TDirection>;

  typedef
    PetscExchange::Handler
    <GridType, VxpeRhsAcquiererAction, TDimension, TDirection>
    VxpeRhsAcquiererHandler;

  using VypeRhsAcquiererAction
          = PetscExchange::VpeRhsAcquiererAction<1, TDimension, TDirection>;
  typedef
    PetscExchange::Handler
    <GridType, VypeRhsAcquiererAction, TDimension, TDirection>
    VypeRhsAcquiererHandler;

  using VzpeRhsAcquiererAction
          = PetscExchange::VpeRhsAcquiererAction<2, TDimension, TDirection>;
  typedef
    PetscExchange::Handler
    <GridType, VzpeRhsAcquiererAction, TDimension, TDirection>
    VzpeRhsAcquiererHandler;
};
}

template <typename T, int D, int U>
IfsfdHandlersBuilder<T, D, U>::
IfsfdHandlersBuilder(Configuration* configuration,
                     SolverType*    simulation)
  : BaseType(configuration, simulation) {}

template <typename T, int D, int U>
void
IfsfdHandlersBuilder<T, D, U>::
setAsInput() {
  this->BaseType::setAsInput();

  typedef
    PetscExchange::VpeInputRhsGenerationAction
    <0, Dimension, Direction>
    VxpeInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VxpeInputRhsGenerationAction, Dimension, Direction>
    VxpeInputRhsGenerationHandler;

  typedef
    PetscExchange::VpeInputRhsGenerationAction
    <1, Dimension, Direction>
    VypeInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VypeInputRhsGenerationAction, Dimension, Direction>
    VypeInputRhsGenerationHandler;

  typedef
    PetscExchange::VpeInputRhsGenerationAction
    <2, Dimension, Direction>
    VzpeInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VzpeInputRhsGenerationAction, Dimension, Direction>
    VzpeInputRhsGenerationHandler;

  using Vpe = Private::VpeTraits<T, D, U>;

  for (int d = 0; d < Dimensions; ++d) {
    if (d == Dimension) {
      if (Direction == 0) {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][Dimension][Direction]
          = Vpe::StencilGenerationHandler::getDirichletLeft(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution());
      } else {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][Dimension][Direction]
          = Vpe::StencilGenerationHandler::getDirichletRight(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution(),
          1);
      }
    } else {
      this->_solver->ghostHandlers()
      ->vpeStencilGeneratorStack[d][Dimension][Direction]
        = Vpe::StencilGenerationHandler::getDirichletMiddle(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution());
    }
  }
  this->_solver->ghostHandlers()
  ->vpeRhsGeneratorStack[0][Dimension][Direction]
    = VxpeInputRhsGenerationHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new VxpeInputRhsGenerationAction(this->_configuration));
  this->_solver->ghostHandlers()
  ->vpeRhsAcquiererStack[0][Dimension][Direction]
    = Vpe::VxpeRhsAcquiererHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new typename Vpe::VxpeRhsAcquiererAction());

  this->_solver->ghostHandlers()
  ->vpeRhsGeneratorStack[1][Dimension][Direction]
    = VypeInputRhsGenerationHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new VypeInputRhsGenerationAction(this->_configuration));
  this->_solver->ghostHandlers()
  ->vpeRhsAcquiererStack[1][Dimension][Direction]
    = Vpe::VypeRhsAcquiererHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new typename Vpe::VypeRhsAcquiererAction());

  if (Dimensions > 2) {
    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[2][Dimension][Direction]
      = VzpeInputRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VzpeInputRhsGenerationAction(this->_configuration));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[2][Dimension][Direction]
      = Vpe::VzpeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new typename Vpe::VzpeRhsAcquiererAction());
  }
}

template <typename T, int D, int U>
void
IfsfdHandlersBuilder<T, D, U>::
setAsParabolicInput() {
  this->BaseType::setAsParabolicInput();

  typedef
    PetscExchange::
    VpeParabolicInputRhsGenerationAction
    <0, Dimension, Direction>
    VxpeParabolicInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VxpeParabolicInputRhsGenerationAction, Dimension, Direction>
    VxpeParabolicInputRhsGenerationHandler;

  typedef
    PetscExchange::
    VpeParabolicInputRhsGenerationAction
    <1, Dimension, Direction>
    VypeParabolicInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VypeParabolicInputRhsGenerationAction, Dimension, Direction>
    VypeParabolicInputRhsGenerationHandler;

  typedef
    PetscExchange::VpeParabolicInputRhsGenerationAction
    <2, Dimension, Direction>
    VzpeParabolicInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VzpeParabolicInputRhsGenerationAction, Dimension, Direction>
    VzpeParabolicInputRhsGenerationHandler;

  using Vpe = Private::VpeTraits<T, D, U>;

  for (int d = 0; d < Dimensions; ++d) {
    if (d == Dimension) {
      if (Direction == 0) {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][Dimension][Direction]
          = Vpe::StencilGenerationHandler::getDirichletLeft(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution());
      } else {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][Dimension][Direction]
          = Vpe::StencilGenerationHandler::getDirichletRight(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution(),
          1);
      }
    } else {
      this->_solver->ghostHandlers()
      ->vpeStencilGeneratorStack[d][Dimension][Direction]
        = Vpe::StencilGenerationHandler::getDirichletMiddle(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution());
    }
  }
  this->_solver->ghostHandlers()
  ->vpeRhsGeneratorStack[0][Dimension][Direction]
    = VxpeParabolicInputRhsGenerationHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new VxpeParabolicInputRhsGenerationAction(this->_configuration));
  this->_solver->ghostHandlers()
  ->vpeRhsAcquiererStack[0][Dimension][Direction]
    = Vpe::VxpeRhsAcquiererHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new typename Vpe::VxpeRhsAcquiererAction());

  this->_solver->ghostHandlers()
  ->vpeRhsGeneratorStack[1][Dimension][Direction]
    = VypeParabolicInputRhsGenerationHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new VypeParabolicInputRhsGenerationAction(this->_configuration));
  this->_solver->ghostHandlers()
  ->vpeRhsAcquiererStack[1][Dimension][Direction]
    = Vpe::VypeRhsAcquiererHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new typename Vpe::VypeRhsAcquiererAction());

  if (Dimensions > 2) {
    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[2][Dimension][Direction]
      = VzpeParabolicInputRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VzpeParabolicInputRhsGenerationAction(this->_configuration));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[2][Dimension][Direction]
      = Vpe::VzpeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new typename Vpe::VzpeRhsAcquiererAction());
  }
}

template <typename T, int D, int U>
void
IfsfdHandlersBuilder<T, D, U>::
setAsOutput() {
  this->BaseType::setAsOutput();

  typedef
    PetscExchange::VpeConstantRhsGenerationAction
    <ScalarType, 0, Dimension, Direction>
    VxpeConstantRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VxpeConstantRhsGenerationAction, Dimension, Direction>
    VxpeConstantRhsGenerationHandler;

  typedef
    PetscExchange::VpeConstantRhsGenerationAction
    <ScalarType, 1, Dimension, Direction>
    VypeConstantRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VypeConstantRhsGenerationAction, Dimension, Direction>
    VypeConstantRhsGenerationHandler;

  typedef
    PetscExchange::VpeConstantRhsGenerationAction
    <ScalarType, 2, Dimension, Direction>
    VzpeConstantRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VzpeConstantRhsGenerationAction, Dimension, Direction>
    VzpeConstantRhsGenerationHandler;

  using Vpe = Private::VpeTraits<T, D, U>;

  for (int d = 0; d < Dimensions; ++d) {
    if (d == Dimension) {
      if (Direction == 0) {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][Dimension][Direction]
          = Vpe::StencilGenerationHandler::getNeumannLeft(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution());
      } else {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][Dimension][Direction]
          = Vpe::StencilGenerationHandler::getNeumannRight(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution(),
          1);
      }
    } else {
      this->_solver->ghostHandlers()
      ->vpeStencilGeneratorStack[d][Dimension][Direction]
        = Vpe::StencilGenerationHandler::getNeumannMiddle(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution());
    }
  }

  this->_solver->ghostHandlers()
  ->vpeRhsGeneratorStack[0][Dimension][Direction]
    = VxpeConstantRhsGenerationHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new VxpeConstantRhsGenerationAction(0.0));
  this->_solver->ghostHandlers()
  ->vpeRhsAcquiererStack[0][Dimension][Direction]
    = Vpe::VxpeRhsAcquiererHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new typename Vpe::VxpeRhsAcquiererAction());

  this->_solver->ghostHandlers()
  ->vpeRhsGeneratorStack[1][Dimension][Direction]
    = VypeConstantRhsGenerationHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new VypeConstantRhsGenerationAction(0.0));
  this->_solver->ghostHandlers()
  ->vpeRhsAcquiererStack[1][Dimension][Direction]
    = Vpe::VypeRhsAcquiererHandler::getHandler(
    this->_solver->memory()->grid(),
    this->_solver->memory()->parallelDistribution(),
    new typename Vpe::VypeRhsAcquiererAction());

  if (Dimensions > 2) {
    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[2][Dimension][Direction]
      = VzpeConstantRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VzpeConstantRhsGenerationAction(0.0));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[2][Dimension][Direction]
      = Vpe::VzpeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new typename Vpe::VzpeRhsAcquiererAction());
  }
}

template <typename T, int D, int U>
void
IfsfdHandlersBuilder<T, D, U>::
setAsMpiExchange() {
  this->BaseType::setAsMpiExchange();
}

template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 0, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 0, 1>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 1, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 1, 1>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 2, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 2>, 2, 1>;

template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 0, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 0, 1>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 1, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 1, 1>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 2, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 2>, 2, 1>;

template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 0, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 0, 1>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 1, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 1, 1>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 2, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<0, double, 3>, 2, 1>;

template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 0, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 0, 1>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 1, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 1, 1>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 2, 0>;
template class IfsfdHandlersBuilder
  <IfsfdSolverTraits<1, double, 3>, 2, 1>;
}
}
}
