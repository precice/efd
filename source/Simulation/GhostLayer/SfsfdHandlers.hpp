#pragma once

#include "InitializationHandler.hpp"
#include "MpiExchangeHandler.hpp"
#include "PetscExchangeHandler.hpp"
#include "PressureStencilHanler.hpp"

namespace Fluid {
namespace Simulation {
namespace GhostLayer {
template <int TDimensions>
class SfsfdHandlers {
public:
  enum {
    Dimensions = TDimensions
  };

  SfsfdHandlers() {
    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        mpiFghExchangeStack[d][d2]
          = MpiExchange::getEmptyFunctor<Dimensions>();
        mpiPressureExchangeStack[d][d2]
          = MpiExchange::getEmptyFunctor<Dimensions>();
        mpiVelocityExchangeStack[d][d2]
          = MpiExchange::getEmptyFunctor<Dimensions>();
        mpiLocationsExchangeStack[d][d2]
          = MpiExchange::getEmptyFunctor<Dimensions>();

        ppeStencilGeneratorStack[d][d2]
          = LsStencilGenerator::getEmptyFunctor<Dimensions>();
        ppeRhsGeneratorStack[d][d2]
          = PetscExchange::getEmptyFunctor<Dimensions>();
        ppeRhsAcquiererStack[d][d2]
          = PetscExchange::getEmptyFunctor<Dimensions>();

        fghInitialization[d][d2]
          = Initialization::getEmptyFunctor<Dimensions>();
        velocityInitialization[d][d2]
          = Initialization::getEmptyFunctor<Dimensions>();
      }
    }
    cornersHandler = [] () {};
  }

  SfsfdHandlers(SfsfdHandlers const& other) = delete;

  ~SfsfdHandlers() {}

  SfsfdHandlers&
  operator=(SfsfdHandlers const& other) = delete;

  void
  executeFghMpiExchange() const {
    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        mpiFghExchangeStack[d][d2]();
      }
    }
  }

  void
  executePressureMpiExchange() const {
    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        mpiPressureExchangeStack[d][d2]();
      }
    }
  }

  void
  executeLocationsMpiExchange() const {
    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        mpiLocationsExchangeStack[d][d2]();
      }
    }
  }

  void
  executeVelocityMpiExchange() const {
    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        mpiVelocityExchangeStack[d][d2]();
      }
    }
    cornersHandler();
  }

  void
  executeFghInitialization() const {
    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        fghInitialization[d][d2]();
      }
    }
  }

  void
  executeVelocityInitialization() const {
    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        velocityInitialization[d][d2]();
      }
    }
  }

  MpiExchange::FunctorStack<Dimensions> mpiFghExchangeStack;
  MpiExchange::FunctorStack<Dimensions> mpiPressureExchangeStack;
  MpiExchange::FunctorStack<Dimensions> mpiVelocityExchangeStack;
  MpiExchange::FunctorStack<Dimensions> mpiLocationsExchangeStack;

  Initialization::FunctorStack<Dimensions> fghInitialization;

  LsStencilGenerator::FunctorStack<Dimensions> ppeStencilGeneratorStack;
  PetscExchange::FunctorStack<Dimensions>      ppeRhsGeneratorStack;
  PetscExchange::FunctorStack<Dimensions>      ppeRhsAcquiererStack;

  Initialization::FunctorStack<Dimensions> velocityInitialization;

  std::function<void()> cornersHandler;
};
}
}
}
