#pragma once

#include "LinearSolver.hpp"
#include "functions.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class SfsfdPeSolver {
public:
  using SolverTraitsType = TSolverTraits;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using GhostHandlersType = typename SolverTraitsType::GhostHandlersType;

  using PpeStencilGeneratorType = PpeStencilGenerator<CellAccessorType,
                                                      ParallelDistributionType>;

  using PpeRhsGeneratorType = PpeRhsGenerator<CellAccessorType>;

  using PpeResultAcquirerType = PpeResultAcquirer1<CellAccessorType>;

  using PpeSolverType
          = LinearSolver<GridType,
                         PpeStencilGeneratorType,
                         PpeRhsGeneratorType,
                         PpeResultAcquirerType>;

  void
  initialize(MemoryType const*        memory,
             GhostHandlersType const* ghost_handlers) {
    _ghostHandlers = ghost_handlers;

    _ppeStencilGenerator.initialize(memory->parallelDistribution());
    _ppeRhsGenerator.initialize(&memory->timeStepSize());

    _ppeSolver.initialize(memory->grid(),
                          memory->parallelDistribution(),
                          &_ppeStencilGenerator,
                          &_ppeRhsGenerator,
                          &_ppeResultAcquierer,
                          &ghost_handlers->ppeStencilGeneratorStack,
                          &ghost_handlers->ppeRhsGeneratorStack,
                          &ghost_handlers->ppeRhsAcquiererStack);
  }

  void
  executeVpe() {
    _ghostHandlers->executeFghInitialization();
  }

  void
  executePpe() {
    _ppeSolver.solve();
  }

private:
  GhostHandlersType const* _ghostHandlers;
  PpeStencilGeneratorType  _ppeStencilGenerator;
  PpeRhsGeneratorType      _ppeRhsGenerator;
  PpeResultAcquirerType    _ppeResultAcquierer;
  PpeSolverType            _ppeSolver;
};

namespace Private {
template <typename TSolverTraits, int TDimensions>
class IfsfdPeSolver3DPart {};
template <typename TSolverTraits>
class IfsfdPeSolver3DPart<TSolverTraits, 2> {
public:
  using SolverTraitsType = TSolverTraits;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using GhostHandlersType = typename SolverTraitsType::GhostHandlersType;

  using VpeStencilGeneratorType
          = VpeStencilGenerator<CellAccessorType,
                                ParallelDistributionType,
                                ParametersType>;

  void
  initialize(MemoryType const*              memory,
             GhostHandlersType const*       ghost_handlers,
             VpeStencilGeneratorType const* stencil_generator) {
    ((void)memory);
    ((void)ghost_handlers);
    ((void)stencil_generator);
  }

  void
  executeVpe() {}
};

template <typename TSolverTraits>
class IfsfdPeSolver3DPart<TSolverTraits, 3> {
public:
  using SolverTraitsType = TSolverTraits;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using GhostHandlersType = typename SolverTraitsType::GhostHandlersType;

  using VpeStencilGeneratorType
          = VpeStencilGenerator<CellAccessorType,
                                ParallelDistributionType,
                                ParametersType>;

  using VzpeRhsGeneratorType = VpeRhsGenerator<CellAccessorType, 2>;

  using VzpeResultAcquirerType = VpeResultAcquirer<2>;

  using VzpeSolverType
          = LinearSolver<GridType,
                         VpeStencilGeneratorType,
                         VzpeRhsGeneratorType,
                         VzpeResultAcquirerType>;

  void
  initialize(MemoryType const*              memory,
             GhostHandlersType const*       ghost_handlers,
             VpeStencilGeneratorType const* stencil_generator) {
    _vzpeSolver.initialize(memory->grid(),
                           memory->parallelDistribution(),
                           stencil_generator,
                           &_vzpeRhsGenerator,
                           &_vzpeResultGenerator,
                           &ghost_handlers->vpeStencilGeneratorStack[2],
                           &ghost_handlers->vpeRhsGeneratorStack[2],
                           &ghost_handlers->vpeRhsAcquiererStack[2]);
  }

  void
  executeVpe() {
    _vzpeSolver.update();
    _vzpeSolver.solve();
  }

private:
  VzpeRhsGeneratorType   _vzpeRhsGenerator;
  VzpeResultAcquirerType _vzpeResultGenerator;
  VzpeSolverType         _vzpeSolver;
};
}

template <typename TSolverTraits>
class IfsfdPeSolver {
public:
  using SolverTraitsType = TSolverTraits;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using GhostHandlersType = typename SolverTraitsType::GhostHandlersType;

  using VpeStencilGeneratorType
          = VpeStencilGenerator<CellAccessorType,
                                ParallelDistributionType,
                                ParametersType>;

  using VxpeRhsGeneratorType = VpeRhsGenerator<CellAccessorType, 0>;

  using VxpeResultAcquirerType = VpeResultAcquirer<0>;

  using VypeRhsGeneratorType = VpeRhsGenerator<CellAccessorType, 1>;

  using VypeResultAcquirerType = VpeResultAcquirer<1>;

  using PpeStencilGeneratorType = PpeStencilGenerator<CellAccessorType,
                                                      ParallelDistributionType>;

  using PpeRhsGeneratorType = PpeRhsGenerator<CellAccessorType>;

  using PpeResultAcquirerType = PpeResultAcquirer2<CellAccessorType>;

  using PpeSolverType
          = LinearSolver<GridType,
                         PpeStencilGeneratorType,
                         PpeRhsGeneratorType,
                         PpeResultAcquirerType>;

  using VxpeSolverType
          = LinearSolver<GridType,
                         VpeStencilGeneratorType,
                         VxpeRhsGeneratorType,
                         VxpeResultAcquirerType>;

  using VypeSolverType
          = LinearSolver<GridType,
                         VpeStencilGeneratorType,
                         VypeRhsGeneratorType,
                         VypeResultAcquirerType>;

  using VzpeSolverType
          = Private::IfsfdPeSolver3DPart<SolverTraitsType, Dimensions>;

  void
  initialize(MemoryType const*        memory,
             GhostHandlersType const* ghost_handlers) {
    _vpeStencilGenerator.initialize(memory->parallelDistribution(),
                                    memory->parameters(),
                                    &memory->timeStepSize());

    _ppeStencilGenerator.initialize(memory->parallelDistribution());

    _ppeRhsGenerator.initialize(&memory->timeStepSize());
    _ppeResultAcquierer.initialize(&memory->timeStepSize());

    _vxpeSolver.initialize(memory->grid(),
                           memory->parallelDistribution(),
                           &_vpeStencilGenerator,
                           &_vxpeRhsGenerator,
                           &_vxpeResultGenerator,
                           &ghost_handlers->vpeStencilGeneratorStack[0],
                           &ghost_handlers->vpeRhsGeneratorStack[0],
                           &ghost_handlers->vpeRhsAcquiererStack[0]);

    _vypeSolver.initialize(memory->grid(),
                           memory->parallelDistribution(),
                           &_vpeStencilGenerator,
                           &_vypeRhsGenerator,
                           &_vypeResultGenerator,
                           &ghost_handlers->vpeStencilGeneratorStack[1],
                           &ghost_handlers->vpeRhsGeneratorStack[1],
                           &ghost_handlers->vpeRhsAcquiererStack[1]);

    _ppeSolver.initialize(memory->grid(),
                          memory->parallelDistribution(),
                          &_ppeStencilGenerator,
                          &_ppeRhsGenerator,
                          &_ppeResultAcquierer,
                          &ghost_handlers->ppeStencilGeneratorStack,
                          &ghost_handlers->ppeRhsGeneratorStack,
                          &ghost_handlers->ppeRhsAcquiererStack);

    _vzpeSolver.initialize(memory,
                           ghost_handlers,
                           &_vpeStencilGenerator);
  }

  void
  executeVpe() {
    _vxpeSolver.update();
    _vxpeSolver.solve();
    _vypeSolver.update();
    _vypeSolver.solve();
    _vzpeSolver.executeVpe();
  }

  void
  executePpe() {
    _ppeSolver.solve();
  }

private:
  VpeStencilGeneratorType _vpeStencilGenerator;

  VxpeRhsGeneratorType   _vxpeRhsGenerator;
  VxpeResultAcquirerType _vxpeResultGenerator;
  VxpeSolverType         _vxpeSolver;

  VypeRhsGeneratorType   _vypeRhsGenerator;
  VypeResultAcquirerType _vypeResultGenerator;
  VypeSolverType         _vypeSolver;

  VzpeSolverType _vzpeSolver;

  PpeStencilGeneratorType _ppeStencilGenerator;
  PpeRhsGeneratorType     _ppeRhsGenerator;
  PpeResultAcquirerType   _ppeResultAcquierer;
  PpeSolverType           _ppeSolver;
};
}
}
