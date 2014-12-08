#ifndef FsiSimulation_GhostCellsHandler_hpp
#define FsiSimulation_GhostCellsHandler_hpp

#include "Solvers/GhostInitializationHandlers.hpp"
#include "Solvers/GhostMpiExchangeHandlers.hpp"
#include "Solvers/GhostPressureStencilHanlers.hpp"
#include "Solvers/GhostPetscExchangeHandlers.hpp"

namespace FsiSimulation {
template <int D>
class GhostCellsHandler {
public:
  GhostCellsHandler() {
    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _pressureStencilStack[d][d2] =
          Solvers::Ghost::PressureStencil::getEmptyFunctor<D>();
        _mpiFghExchangeStack[d][d2] =
          Solvers::Ghost::MpiExchange::getEmptyFunctor<D>();
        _mpiPressureExchangeStack[d][d2] =
          Solvers::Ghost::MpiExchange::getEmptyFunctor<D>();
        _mpiVelocityExchangeStack[d][d2] =
          Solvers::Ghost::MpiExchange::getEmptyFunctor<D>();
        _fghInitialization[d][d2] =
          Solvers::Ghost::Initialization::getEmptyFunctor<D>();
        _rhsInitialization[d][d2] =
          Solvers::Ghost::PetscExchange::getEmptyFunctor<D>();
        _pressureInitialization[d][d2] =
          Solvers::Ghost::PetscExchange::getEmptyFunctor<D>();
        _velocityInitialization[d][d2] =
          Solvers::Ghost::Initialization::getEmptyFunctor<D>();
      }
    }
  }

  GhostCellsHandler(GhostCellsHandler const& other) = delete;

  ~GhostCellsHandler() {}

  GhostCellsHandler&
  operator=(GhostCellsHandler const& other) = delete;

  Solvers::Ghost::PressureStencil::FunctorStack<D>   _pressureStencilStack;
  Solvers::Ghost::MpiExchange::FunctorStack<D>       _mpiFghExchangeStack;
  Solvers::Ghost::MpiExchange::FunctorStack<D>       _mpiVelocityExchangeStack;
  Solvers::Ghost::MpiExchange::FunctorStack<D>       _mpiPressureExchangeStack;
  Solvers::Ghost::Initialization::FunctorStack<D>    _fghInitialization;
  Solvers::Ghost::PetscExchange::FunctorStack<D> _rhsInitialization;
  Solvers::Ghost::PetscExchange::FunctorStack<D> _pressureInitialization;
  Solvers::Ghost::Initialization::FunctorStack<D>    _velocityInitialization;

private:
};
}
#endif
