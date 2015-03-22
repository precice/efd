#pragma once

#include "GridGeometry.hpp"
#include "SolverTraits.hpp"

#include "GhostLayer/SfsfdHandlersBuilder.hpp"
#include "GhostLayer/IfsfdHandlersBuilder.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TScalar,
          int TDimensions,
          int TSolverType>
struct SolverBuilderTraits {};

template <typename TScalar, int TDimensions>
struct SolverBuilderTraits<TScalar, TDimensions, 0> {
  using SolverTraitsType = SfsfdSolverTraits
                           <UniformGridGeometry<TScalar, TDimensions>,
                            TScalar,
                            TDimensions>;

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = GhostLayer::SfsfdHandlersBuilder
            <SolverTraitsType, TDimension, TDirection>;
};

template <typename TScalar, int TDimensions>
struct SolverBuilderTraits<TScalar, TDimensions, 1> {
  using SolverTraitsType = IfsfdSolverTraits
                           <UniformGridGeometry<TScalar, TDimensions>,
                            TScalar,
                            TDimensions>;

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = GhostLayer::IfsfdHandlersBuilder
            <SolverTraitsType, TDimension, TDirection>;
};
}
}
