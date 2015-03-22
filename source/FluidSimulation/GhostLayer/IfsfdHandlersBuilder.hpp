#pragma once

#include "FsfdHandlersBuilder.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <typename TSolverTraits,
          int TDimension,
          int TDirection>
struct IfsfdHandlersBuilderTraits {
  using SolverTraitsType = TSolverTraits;
  enum {
    Dimension  = TDimension,
    Direction  = TDirection,
    Dimensions = SolverTraitsType::Dimensions
  };

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  static ScalarType*
  projectionTermAccessor(CellAccessorType const& accessor) {
    return &accessor.projectionTerm();
  }

  using
    PressureMpiExchangeHandler
      = MpiExchange::Handler
        <ScalarType,
         typename GridType::BaseType,
         IfsfdHandlersBuilderTraits::projectionTermAccessor,
         TDimension,
         TDirection>;

  using PpeRhsAcquiererAction
          = PetscExchange::PpeRhsAcquiererAction2;
};

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class IfsfdHandlersBuilder
  : public FsfdHandlersBuilder
    < IfsfdHandlersBuilderTraits < TSolverTraits, TDimension, TDirection >> {
public:
  using HandlersBuilderTraitsType
          = IfsfdHandlersBuilderTraits<TSolverTraits, TDimension, TDirection>;

  using BaseType = FsfdHandlersBuilder<HandlersBuilderTraitsType>;

  using SolverTraitsType = TSolverTraits;

  using SolverType = typename SolverTraitsType::SolverType;

  enum {
    Dimension  = TDimension,
    Direction  = TDirection,
    Dimensions = SolverTraitsType::Dimensions
  };

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  typedef
    typename HandlersBuilderTraitsType::PressureMpiExchangeHandler
    PressureMpiExchangeHandler;

  typedef
    typename HandlersBuilderTraitsType::PpeRhsAcquiererAction
    PpeRhsAcquiererAction;

  typedef
    LsStencilGenerator::Handler
    <GridType, TDimension, TDirection>
    VpeStencilGenerationHandler;

  typedef
    PetscExchange::ConstantRhsGenerationAction
    <CellAccessorType, TDimension>
    VpeConstantRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VpeConstantRhsGenerationAction, TDimension, TDirection>
    VpeConstantRhsGenerationHandler;

  typedef
    PetscExchange::VpeInputRhsGenerationAction
    <0, TDimension, TDirection>
    VxpeInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VxpeInputRhsGenerationAction, TDimension, TDirection>
    VxpeInputRhsGenerationHandler;

  typedef
    PetscExchange::
    VpeParabolicInputRhsGenerationAction
    <0, TDimension, TDirection>
    VxpeParabolicInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VxpeParabolicInputRhsGenerationAction, TDimension, TDirection>
    VxpeParabolicInputRhsGenerationHandler;

  typedef
    PetscExchange::VpeRhsAcquiererAction<0>
    VxpeRhsAcquiererAction;
  typedef
    PetscExchange::Handler
    <GridType, VxpeRhsAcquiererAction, TDimension, TDirection>
    VxpeRhsAcquiererHandler;

  typedef
    PetscExchange::VpeInputRhsGenerationAction
    <1, TDimension, TDirection>
    VypeInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VypeInputRhsGenerationAction, TDimension, TDirection>
    VypeInputRhsGenerationHandler;

  typedef
    PetscExchange::
    VpeParabolicInputRhsGenerationAction
    <1, TDimension, TDirection>
    VypeParabolicInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VypeParabolicInputRhsGenerationAction, TDimension, TDirection>
    VypeParabolicInputRhsGenerationHandler;

  typedef
    PetscExchange::VpeRhsAcquiererAction<1>
    VypeRhsAcquiererAction;
  typedef
    PetscExchange::Handler
    <GridType, VypeRhsAcquiererAction, TDimension, TDirection>
    VypeRhsAcquiererHandler;

  typedef
    PetscExchange::VpeInputRhsGenerationAction
    <2, TDimension, TDirection>
    VzpeInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VzpeInputRhsGenerationAction, TDimension, TDirection>
    VzpeInputRhsGenerationHandler;

  typedef
    PetscExchange::VpeParabolicInputRhsGenerationAction
    <2, TDimension, TDirection>
    VzpeParabolicInputRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VzpeParabolicInputRhsGenerationAction, TDimension, TDirection>
    VzpeParabolicInputRhsGenerationHandler;

  typedef
    PetscExchange::VpeRhsAcquiererAction<2>
    VzpeRhsAcquiererAction;
  typedef
    PetscExchange::Handler
    <GridType, VzpeRhsAcquiererAction, TDimension, TDirection>
    VzpeRhsAcquiererHandler;

public:
  IfsfdHandlersBuilder(Configuration* configuration,
                       SolverType*    simulation)
    : BaseType(configuration, simulation) {}

  inline void
  setAsInput() {
    this->BaseType::setAsInput();

    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          this->_solver->ghostHandlers()
          ->vpeStencilGeneratorStack[d][TDimension][TDirection]
            = VpeStencilGenerationHandler::getDirichletLeft(
            this->_solver->memory()->grid(),
            this->_solver->memory()->parallelDistribution());
        } else {
          this->_solver->ghostHandlers()
          ->vpeStencilGeneratorStack[d][TDimension][TDirection]
            = VpeStencilGenerationHandler::getDirichletRight(
            this->_solver->memory()->grid(),
            this->_solver->memory()->parallelDistribution());
        }
      } else {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][TDimension][TDirection]
          = VpeStencilGenerationHandler::getDirichletMiddle(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution());
      }
    }
    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[0][TDimension][TDirection]
      = VxpeInputRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VxpeInputRhsGenerationAction(this->_configuration));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[0][TDimension][TDirection]
      = VxpeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VxpeRhsAcquiererAction());

    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[1][TDimension][TDirection]
      = VypeInputRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VypeInputRhsGenerationAction(this->_configuration));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[1][TDimension][TDirection] =
      VypeRhsAcquiererHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VypeRhsAcquiererAction());

    if (Dimensions > 2) {
      this->_solver->ghostHandlers()
      ->vpeRhsGeneratorStack[2][TDimension][TDirection]
        = VzpeInputRhsGenerationHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VzpeInputRhsGenerationAction(this->_configuration));
      this->_solver->ghostHandlers()
      ->vpeRhsAcquiererStack[2][TDimension][TDirection]
        = VzpeRhsAcquiererHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VzpeRhsAcquiererAction());
    }
  }

  inline void
  setAsParabolicInput() {
    this->BaseType::setAsParabolicInput();

    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          this->_solver->ghostHandlers()
          ->vpeStencilGeneratorStack[d][TDimension][TDirection]
            = VpeStencilGenerationHandler::getDirichletLeft(
            this->_solver->memory()->grid(),
            this->_solver->memory()->parallelDistribution());
        } else {
          this->_solver->ghostHandlers()
          ->vpeStencilGeneratorStack[d][TDimension][TDirection]
            = VpeStencilGenerationHandler::getDirichletRight(
            this->_solver->memory()->grid(),
            this->_solver->memory()->parallelDistribution());
        }
      } else {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][TDimension][TDirection]
          = VpeStencilGenerationHandler::getDirichletMiddle(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution());
      }
    }
    this->_solver->ghostHandlers()->vpeRhsGeneratorStack[0][TDimension][
      TDirection] =
      VxpeParabolicInputRhsGenerationHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VxpeParabolicInputRhsGenerationAction(this->_configuration));
    this->_solver->ghostHandlers()->vpeRhsAcquiererStack[0][TDimension][
      TDirection] =
      VxpeRhsAcquiererHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VxpeRhsAcquiererAction());

    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[1][TDimension][TDirection]
      = VypeParabolicInputRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VypeParabolicInputRhsGenerationAction(this->_configuration));
    this->_solver->ghostHandlers()->vpeRhsAcquiererStack[1][TDimension][
      TDirection] =
      VypeRhsAcquiererHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VypeRhsAcquiererAction());

    if (Dimensions > 2) {
      this->_solver->ghostHandlers()
      ->vpeRhsGeneratorStack[2][TDimension][TDirection] =
        VzpeParabolicInputRhsGenerationHandler::getHandler(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution(),
          new VzpeParabolicInputRhsGenerationAction(this->_configuration));
      this->_solver->ghostHandlers()
      ->vpeRhsAcquiererStack[2][TDimension][TDirection]
        = VzpeRhsAcquiererHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VzpeRhsAcquiererAction());
    }
  }

  inline void
  setAsOutput() {
    this->BaseType::setAsOutput();

    for (int d = 0; d < Dimensions; ++d) {
      if (d == TDimension) {
        if (TDirection == 0) {
          this->_solver->ghostHandlers()
          ->vpeStencilGeneratorStack[d][TDimension][TDirection]
            = VpeStencilGenerationHandler::getNeumannLeft(
            this->_solver->memory()->grid(),
            this->_solver->memory()->parallelDistribution());
        } else {
          this->_solver->ghostHandlers()
          ->vpeStencilGeneratorStack[d][TDimension][TDirection]
            = VpeStencilGenerationHandler::getNeumannRight(
            this->_solver->memory()->grid(),
            this->_solver->memory()->parallelDistribution());
        }
      } else {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][TDimension][TDirection]
          = VpeStencilGenerationHandler::getNeumannMiddle(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution());
      }
    }

    int offset = 0;

    if (TDimension == 0) {
      if (TDirection == 1) {
        offset = 1;
      }
    }
    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[0][TDimension][TDirection]
      = VpeConstantRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VpeConstantRhsGenerationAction(0.0, offset));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[0][TDimension][TDirection]
      = VxpeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VxpeRhsAcquiererAction());

    offset = 0;

    if (TDimension == 1) {
      if (TDirection == 1) {
        offset = 1;
      }
    }

    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[1][TDimension][TDirection]
      = VpeConstantRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VpeConstantRhsGenerationAction(0.0, offset));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[1][TDimension][TDirection]
      = VypeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VypeRhsAcquiererAction());

    if (Dimensions > 2) {
      offset = 0;

      if (TDimension == 2) {
        if (TDirection == 1) {
          offset = 1;
        }
      }
      this->_solver->ghostHandlers()
      ->vpeRhsGeneratorStack[2][TDimension][TDirection]
        = VpeConstantRhsGenerationHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VpeConstantRhsGenerationAction(0.0, offset));
      this->_solver->ghostHandlers()
      ->vpeRhsAcquiererStack[2][TDimension][TDirection]
        = VzpeRhsAcquiererHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VzpeRhsAcquiererAction());
    }
  }

  void
  setAsMpiExchange() {
    this->BaseType::setAsMpiExchange();
  }
};
}
}
}
