#ifndef FsiSimulation_GhostCellsHandler_hpp
#define FsiSimulation_GhostCellsHandler_hpp

#include "InitializationHandler.hpp"
#include "MpiExchangeHandler.hpp"
#include "PressureStencilHanler.hpp"
#include "PetscExchangeHandler.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template<int D>
class Handlers {
public:
  Handlers() {
    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _pressureStencilStack[d][d2] =
        PressureStencil::getEmptyFunctor<D>();
        _mpiFghExchangeStack[d][d2] =
        MpiExchange::getEmptyFunctor<D>();
        _mpiPressureExchangeStack[d][d2] =
        MpiExchange::getEmptyFunctor<D>();
        _mpiVelocityExchangeStack[d][d2] =
        MpiExchange::getEmptyFunctor<D>();
        _fghInitialization[d][d2] =
       Initialization::getEmptyFunctor<D>();
        _rhsInitialization[d][d2] =
        PetscExchange::getEmptyFunctor<D>();
        _pressureInitialization[d][d2] =
      PetscExchange::getEmptyFunctor<D>();
        _velocityInitialization[d][d2] =
       Initialization::getEmptyFunctor<D>();
      }
    }
  }

  Handlers(Handlers const& other) = delete;

  ~Handlers() {
  }

  Handlers&
  operator=(Handlers const& other) = delete;

  PressureStencil::FunctorStack<D> _pressureStencilStack;
  MpiExchange::FunctorStack<D> _mpiFghExchangeStack;
  MpiExchange::FunctorStack<D> _mpiVelocityExchangeStack;
  MpiExchange::FunctorStack<D> _mpiPressureExchangeStack;
  Initialization::FunctorStack<D> _fghInitialization;
  PetscExchange::FunctorStack<D> _rhsInitialization;
  PetscExchange::FunctorStack<D> _pressureInitialization;
  Initialization::FunctorStack<D> _velocityInitialization;

private:
};
}
}
}
#endif
