#pragma once

#include "SfsfdHandlers.hpp"

namespace Fluid {
namespace Simulation {
namespace GhostLayer {
template <int TDimensions>
class IfsfdHandlers : public SfsfdHandlers<TDimensions> {
public:
  using BaseType = SfsfdHandlers<TDimensions>;

  enum {
    Dimensions = BaseType::Dimensions
  };

  IfsfdHandlers() : BaseType() {
    for (int i = 0; i < Dimensions; ++i) {
      for (int d = 0; d < Dimensions; ++d) {
        for (int d2 = 0; d2 < 2; ++d2) {
          vpeStencilGeneratorStack[i][d][d2] =
            LsStencilGenerator::getEmptyFunctor<Dimensions>();
          vpeRhsGeneratorStack[i][d][d2] =
            PetscExchange::getEmptyFunctor<Dimensions>();
          vpeRhsAcquiererStack[i][d][d2] =
            PetscExchange::getEmptyFunctor<Dimensions>();
        }
      }
    }
  }

  IfsfdHandlers(IfsfdHandlers const& other) = delete;

  ~IfsfdHandlers() {}

  IfsfdHandlers&
  operator=(IfsfdHandlers const& other) = delete;

  std::array<LsStencilGenerator::FunctorStack<Dimensions>, Dimensions>
  vpeStencilGeneratorStack;

  std::array<PetscExchange::FunctorStack<Dimensions>, Dimensions>
  vpeRhsGeneratorStack;

  std::array<PetscExchange::FunctorStack<Dimensions>, Dimensions>
  vpeRhsAcquiererStack;
};
}
}
}
