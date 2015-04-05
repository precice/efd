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
  getProjectionTerm(CellAccessorType const& accessor) {
    return &accessor.projectionTerm();
  }

  static void
  setProjectionTerm(CellAccessorType const& accessor,
                    int const& index,
                    ScalarType const& value) {
    ((void)index);
    accessor.pressure() += value;
    accessor.projectionTerm() = value;
  }

  using PressureMpiExchangeHandler
      = MpiExchange::Handler
        <ScalarType,
         1,
         typename GridType::BaseType,
         IfsfdHandlersBuilderTraits::getProjectionTerm,
         IfsfdHandlersBuilderTraits::setProjectionTerm,
         TDimension,
         TDirection>;

  using PpeRhsAcquiererAction
          = PetscExchange::PpeRhsAcquiererAction2;
};

template <typename TSolverTraits,
          int TDimension,
          int TDirection>
class IfsfdHandlersBuilder :
  public FsfdHandlersBuilder
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
    PetscExchange::VpeConstantRhsGenerationAction
    <ScalarType, 0, TDimension, TDirection>
    VxpeConstantRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VxpeConstantRhsGenerationAction, TDimension, TDirection>
    VxpeConstantRhsGenerationHandler;

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

  using VxpeRhsAcquiererAction
          = PetscExchange::VpeRhsAcquiererAction<0, TDimension, TDirection>;

  typedef
    PetscExchange::Handler
    <GridType, VxpeRhsAcquiererAction, TDimension, TDirection>
    VxpeRhsAcquiererHandler;

  typedef
    PetscExchange::VpeConstantRhsGenerationAction
    <ScalarType, 1, TDimension, TDirection>
    VypeConstantRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VypeConstantRhsGenerationAction, TDimension, TDirection>
    VypeConstantRhsGenerationHandler;

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

  using VypeRhsAcquiererAction
          = PetscExchange::VpeRhsAcquiererAction<1, TDimension, TDirection>;
  typedef
    PetscExchange::Handler
    <GridType, VypeRhsAcquiererAction, TDimension, TDirection>
    VypeRhsAcquiererHandler;

  typedef
    PetscExchange::VpeConstantRhsGenerationAction
    <ScalarType, 2, TDimension, TDirection>
    VzpeConstantRhsGenerationAction;
  typedef
    PetscExchange::Handler
    <GridType, VzpeConstantRhsGenerationAction, TDimension, TDirection>
    VzpeConstantRhsGenerationHandler;

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

  using VzpeRhsAcquiererAction
          = PetscExchange::VpeRhsAcquiererAction<2, TDimension, TDirection>;
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
            this->_solver->memory()->parallelDistribution(),
            1);
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
    ->vpeRhsAcquiererStack[1][TDimension][TDirection]
      = VypeRhsAcquiererHandler::getHandler(
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
            this->_solver->memory()->parallelDistribution(),
            1);
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
      = VxpeParabolicInputRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VxpeParabolicInputRhsGenerationAction(this->_configuration));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[0][TDimension][TDirection]
      = VxpeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VxpeRhsAcquiererAction());

    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[1][TDimension][TDirection]
      = VypeParabolicInputRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VypeParabolicInputRhsGenerationAction(this->_configuration));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[1][TDimension][TDirection]
      = VypeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VypeRhsAcquiererAction());

    if (Dimensions > 2) {
      this->_solver->ghostHandlers()
      ->vpeRhsGeneratorStack[2][TDimension][TDirection]
        = VzpeParabolicInputRhsGenerationHandler::getHandler(
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
            this->_solver->memory()->parallelDistribution(),
            1);
        }
      } else {
        this->_solver->ghostHandlers()
        ->vpeStencilGeneratorStack[d][TDimension][TDirection]
          = VpeStencilGenerationHandler::getNeumannMiddle(
          this->_solver->memory()->grid(),
          this->_solver->memory()->parallelDistribution());
      }
    }

    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[0][TDimension][TDirection]
      = VxpeConstantRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VxpeConstantRhsGenerationAction(0.0));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[0][TDimension][TDirection]
      = VxpeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VxpeRhsAcquiererAction());

    this->_solver->ghostHandlers()
    ->vpeRhsGeneratorStack[1][TDimension][TDirection]
      = VypeConstantRhsGenerationHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VypeConstantRhsGenerationAction(0.0));
    this->_solver->ghostHandlers()
    ->vpeRhsAcquiererStack[1][TDimension][TDirection]
      = VypeRhsAcquiererHandler::getHandler(
      this->_solver->memory()->grid(),
      this->_solver->memory()->parallelDistribution(),
      new VypeRhsAcquiererAction());

    if (Dimensions > 2) {
      this->_solver->ghostHandlers()
      ->vpeRhsGeneratorStack[2][TDimension][TDirection]
        = VzpeConstantRhsGenerationHandler::getHandler(
        this->_solver->memory()->grid(),
        this->_solver->memory()->parallelDistribution(),
        new VzpeConstantRhsGenerationAction(0.0));
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
