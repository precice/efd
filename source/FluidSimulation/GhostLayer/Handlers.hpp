#ifndef FsiSimulation_FluidSimulation_GhostLayer_Handler_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_Handler_hpp

#include "InitializationHandler.hpp"
#include "MpiExchangeHandler.hpp"
#include "PetscExchangeHandler.hpp"
#include "PressureStencilHanler.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <int TDimensions>
class Handlers {
public:
 enum {
   Dimensions = TDimensions
 };

  Handlers() {
    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        mpiFghExchangeStack[d][d2] =
          MpiExchange::getEmptyFunctor<Dimensions>();
        mpiPressureExchangeStack[d][d2] =
          MpiExchange::getEmptyFunctor<Dimensions>();
        mpiVelocityExchangeStack[d][d2] =
          MpiExchange::getEmptyFunctor<Dimensions>();
        fghInitialization[d][d2] =
          Initialization::getEmptyFunctor<Dimensions>();

        ppeStencilGeneratorStack[d][d2] =
          LsStencilGenerator::getEmptyFunctor<Dimensions>();
        ppeRhsGeneratorStack[d][d2] =
          PetscExchange::getEmptyFunctor<Dimensions>();
        ppeRhsAcquiererStack[d][d2] =
          PetscExchange::getEmptyFunctor<Dimensions>();

        for (int i = 0; i < Dimensions; ++i) {
          vpeStencilGeneratorStack[i][d][d2] =
            LsStencilGenerator::getEmptyFunctor<Dimensions>();
          vpeRhsGeneratorStack[i][d][d2] =
            PetscExchange::getEmptyFunctor<Dimensions>();
          vpeRhsAcquiererStack[i][d][d2] =
            PetscExchange::getEmptyFunctor<Dimensions>();
        }

        velocityInitialization[d][d2] =
          Initialization::getEmptyFunctor<Dimensions>();
      }
    }
  }

  Handlers(Handlers const& other) = delete;

  ~Handlers() {}

  Handlers&
  operator=(Handlers const& other) = delete;

  MpiExchange::FunctorStack<Dimensions> mpiFghExchangeStack;
  MpiExchange::FunctorStack<Dimensions> mpiVelocityExchangeStack;
  MpiExchange::FunctorStack<Dimensions> mpiPressureExchangeStack;

  Initialization::FunctorStack<Dimensions> fghInitialization;

  LsStencilGenerator::FunctorStack<Dimensions> ppeStencilGeneratorStack;
  PetscExchange::FunctorStack<Dimensions>      ppeRhsGeneratorStack;
  PetscExchange::FunctorStack<Dimensions>      ppeRhsAcquiererStack;

  std::array<LsStencilGenerator::FunctorStack<Dimensions>, Dimensions> vpeStencilGeneratorStack;
  std::array<PetscExchange::FunctorStack<Dimensions>, Dimensions>      vpeRhsGeneratorStack;
  std::array<PetscExchange::FunctorStack<Dimensions>, Dimensions>      vpeRhsAcquiererStack;

  Initialization::FunctorStack<Dimensions> velocityInitialization;

private:
};
}
}
}
#endif
