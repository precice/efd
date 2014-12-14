#ifndef FsiSimulation_FluidSimulation_GhostLayer_Handler_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_Handler_hpp

#include "InitializationHandler.hpp"
#include "MpiExchangeHandler.hpp"
#include "PetscExchangeHandler.hpp"
#include "PressureStencilHanler.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <int TD>
class Handlers {
public:
  Handlers() {
    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        pressureStencilStack[d][d2] =
          PressureStencil::getEmptyFunctor<TD>();
        mpiFghExchangeStack[d][d2] =
          MpiExchange::getEmptyFunctor<TD>();
        mpiPressureExchangeStack[d][d2] =
          MpiExchange::getEmptyFunctor<TD>();
        mpiVelocityExchangeStack[d][d2] =
          MpiExchange::getEmptyFunctor<TD>();
        fghInitialization[d][d2] =
          Initialization::getEmptyFunctor<TD>();
        rhsInitialization[d][d2] =
          PetscExchange::getEmptyFunctor<TD>();
        pressureInitialization[d][d2] =
          PetscExchange::getEmptyFunctor<TD>();
        velocityInitialization[d][d2] =
          Initialization::getEmptyFunctor<TD>();
      }
    }
  }

  Handlers(Handlers const& other) = delete;

  ~Handlers() {}

  Handlers&
  operator=(Handlers const& other) = delete;

  PressureStencil::FunctorStack<TD> pressureStencilStack;
  MpiExchange::FunctorStack<TD>     mpiFghExchangeStack;
  MpiExchange::FunctorStack<TD>     mpiVelocityExchangeStack;
  MpiExchange::FunctorStack<TD>     mpiPressureExchangeStack;
  Initialization::FunctorStack<TD>  fghInitialization;
  PetscExchange::FunctorStack<TD>   rhsInitialization;
  PetscExchange::FunctorStack<TD>   pressureInitialization;
  Initialization::FunctorStack<TD>  velocityInitialization;

private:
};
}
}
}
#endif
