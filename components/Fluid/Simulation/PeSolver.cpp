#include "PeSolver.hpp"

#include "GhostLayer/IfsfdHandlers.hpp"
#include "GhostLayer/SfsfdHandlers.hpp"
#include "ParallelDistribution.hpp"
#include "Grid.hpp"
#include "GridGeometry.hpp"
#include "IfsfdCellAccessor.hpp"
#include "IfsfdMemory.hpp"
#include "PetscLinearSolver.hpp"
#include "PetscLinearSolverDefinitions.hpp"
#include "SfsfdCellAccessor.hpp"
#include "SfsfdMemory.hpp"
#include "SolverTraits.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename T>
void
SfsfdPeSolver<T>::
initialize(MemoryType const*        memory,
           GhostHandlersType const* ghost_handlers) {
  _ghostHandlers = ghost_handlers;

  using PpeStencilGeneratorType = PpeStencilGenerator<SolverTraitsType>;

  using PpeRhsGeneratorType = PpeRhsGenerator<SolverTraitsType>;

  using PpeResultAcquirerType = PpeResultAcquirer1;

  using PpeSolverType = PetscLinearSolver<GridType,
                                          PpeStencilGeneratorType,
                                          PpeRhsGeneratorType,
                                          PpeResultAcquirerType>;

  std::unique_ptr<PpeSolverType> ppeSolver(new PpeSolverType());

  ppeSolver->stencilGenerator.initialize(memory->parallelDistribution());
  ppeSolver->rhsGenerator.initialize(&memory->timeStepSize());

  ppeSolver->initialize(memory->grid(),
                        memory->parallelDistribution(),
                        &ghost_handlers->ppeStencilGeneratorStack,
                        &ghost_handlers->ppeRhsGeneratorStack,
                        &ghost_handlers->ppeRhsAcquiererStack);
  _ppeSolver.reset(ppeSolver.release());
  _ppeSolver->update();
}

template <typename T>
void
SfsfdPeSolver<T>::
executeVpe() {
  _ghostHandlers->executeFghInitialization();
}

template <typename T>
void
SfsfdPeSolver<T>::
executePpe() {
  _ppeSolver->solve();
}

namespace Private {
template <typename T>
void
IfsfdPeSolver3DPart<T, 3>::
initialize(MemoryType const*        memory,
           GhostHandlersType const* ghost_handlers) {
  _memory = memory;

  using VpeStencilGeneratorType = VpeStencilGenerator<SolverTraitsType>;

  using VzpeRhsGeneratorType = VpeRhsGenerator<2>;

  using VzpeResultAcquirerType = VpeResultAcquirer<2>;

  using VzpeSolverType
          = PetscLinearSolver<GridType,
                              VpeStencilGeneratorType,
                              VzpeRhsGeneratorType,
                              VzpeResultAcquirerType>;

  std::unique_ptr<VzpeSolverType> vzpeSolver(new VzpeSolverType());

  vzpeSolver->initialize(memory->grid(),
                         memory->parallelDistribution(),
                         &ghost_handlers->vpeStencilGeneratorStack[2],
                         &ghost_handlers->vpeRhsGeneratorStack[2],
                         &ghost_handlers->vpeRhsAcquiererStack[2]);
  _vzpeSolver.reset(vzpeSolver.release());

  _vzpeAbsoluteTolerance = _vzpeSolver->absoluteTolerance();
  _vzpeRelativeTolerance = _vzpeSolver->relativeTolerance();
}

template <typename T>
void
IfsfdPeSolver3DPart<T, 3>::
executeVpe() {
  _vzpeSolver->absoluteTolerance(
    _memory->timeStepSize() * _vzpeAbsoluteTolerance / 1000.0);
  _vzpeSolver->relativeTolerance(
    _memory->timeStepSize() * _vzpeRelativeTolerance / 1000.0);
  _vzpeSolver->update();
  _vzpeSolver->solve();
}
}

template <typename T>
void
IfsfdPeSolver<T>::
initialize(MemoryType const*        memory,
           GhostHandlersType const* ghost_handlers) {
  _memory = memory;

  using VpeStencilGeneratorType = VpeStencilGenerator<SolverTraitsType>;

  using VxpeRhsGeneratorType = VpeRhsGenerator<0>;

  using VxpeResultAcquirerType = VpeResultAcquirer<0>;

  using VypeRhsGeneratorType = VpeRhsGenerator<1>;

  using VypeResultAcquirerType = VpeResultAcquirer<1>;

  using PpeStencilGeneratorType = PpeStencilGenerator<SolverTraitsType>;

  using PpeRhsGeneratorType = PpeRhsGenerator<SolverTraitsType>;

  using PpeResultAcquirerType = PpeResultAcquirer2<SolverTraitsType>;

  using PpeSolverType
          = PetscLinearSolver<GridType,
                              PpeStencilGeneratorType,
                              PpeRhsGeneratorType,
                              PpeResultAcquirerType>;

  using VxpeSolverType
          = PetscLinearSolver<GridType,
                              VpeStencilGeneratorType,
                              VxpeRhsGeneratorType,
                              VxpeResultAcquirerType>;

  using VypeSolverType
          = PetscLinearSolver<GridType,
                              VpeStencilGeneratorType,
                              VypeRhsGeneratorType,
                              VypeResultAcquirerType>;

  std::unique_ptr<PpeSolverType>  ppeSolver(new PpeSolverType());
  std::unique_ptr<VxpeSolverType> vxpeSolver(new VxpeSolverType());
  std::unique_ptr<VypeSolverType> vypeSolver(new VypeSolverType());

  ppeSolver->stencilGenerator.initialize(memory->parallelDistribution());

  ppeSolver->rhsGenerator.initialize(&memory->timeStepSize());
  ppeSolver->resultAcquierer.initialize(&memory->timeStepSize());

  vxpeSolver->stencilGenerator.initialize(memory->parallelDistribution(),
                                          memory->parameters(),
                                          &memory->timeStepSize());

  vypeSolver->stencilGenerator.initialize(memory->parallelDistribution(),
                                          memory->parameters(),
                                          &memory->timeStepSize());

  ppeSolver->initialize(memory->grid(),
                        memory->parallelDistribution(),
                        &ghost_handlers->ppeStencilGeneratorStack,
                        &ghost_handlers->ppeRhsGeneratorStack,
                        &ghost_handlers->ppeRhsAcquiererStack);
  _ppeSolver.reset(ppeSolver.release());
  _ppeSolver->update();

  vxpeSolver->initialize(memory->grid(),
                         memory->parallelDistribution(),
                         &ghost_handlers->vpeStencilGeneratorStack[0],
                         &ghost_handlers->vpeRhsGeneratorStack[0],
                         &ghost_handlers->vpeRhsAcquiererStack[0]);
  _vxpeSolver.reset(vxpeSolver.release());

  vypeSolver->initialize(memory->grid(),
                         memory->parallelDistribution(),
                         &ghost_handlers->vpeStencilGeneratorStack[1],
                         &ghost_handlers->vpeRhsGeneratorStack[1],
                         &ghost_handlers->vpeRhsAcquiererStack[1]);
  _vypeSolver.reset(vypeSolver.release());

  _vzpeSolver.initialize(memory,
                         ghost_handlers);

  _vxpeAbsoluteTolerance = _vxpeSolver->absoluteTolerance();
  _vxpeRelativeTolerance = _vxpeSolver->relativeTolerance();
  _vypeAbsoluteTolerance = _vypeSolver->absoluteTolerance();
  _vypeRelativeTolerance = _vypeSolver->relativeTolerance();
}

template <typename T>
void
IfsfdPeSolver<T>::
executeVpe() {
  _vxpeSolver->absoluteTolerance(
    _memory->timeStepSize() * _vxpeAbsoluteTolerance / 1000.0);
  _vxpeSolver->relativeTolerance(
    _memory->timeStepSize() * _vxpeRelativeTolerance / 1000.0);
  _vxpeSolver->update();
  _vxpeSolver->solve();

  _vypeSolver->absoluteTolerance(
    _memory->timeStepSize() * _vypeAbsoluteTolerance / 1000.0);
  _vypeSolver->relativeTolerance(
    _memory->timeStepSize() * _vypeRelativeTolerance / 1000.0);
  _vypeSolver->update();
  _vypeSolver->solve();
  _vzpeSolver.executeVpe();
}

template <typename T>
void
IfsfdPeSolver<T>::
executePpe() {
  _ppeSolver->solve();
}

template class SfsfdPeSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 0, double, 2 >>;
template class SfsfdPeSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 1, double, 2 >>;
template class SfsfdPeSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 0, double, 2 >>;
template class SfsfdPeSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 1, double, 2 >>;
template class SfsfdPeSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 0, double, 3 >>;
template class SfsfdPeSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 1, double, 3 >>;
template class SfsfdPeSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 0, double, 3 >>;
template class SfsfdPeSolver
  < SfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 1, double, 3 >>;

template class IfsfdPeSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 0, double, 2 >>;
template class IfsfdPeSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 0, 1, double, 2 >>;
template class IfsfdPeSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 0, double, 2 >>;
template class IfsfdPeSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 2>, 1, 1, double, 2 >>;
template class IfsfdPeSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 0, double, 3 >>;
template class IfsfdPeSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 0, 1, double, 3 >>;
template class IfsfdPeSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 0, double, 3 >>;
template class IfsfdPeSolver
  < IfsfdSolverTraits < UniformGridGeometry<double, 3>, 1, 1, double, 3 >>;
}
}
