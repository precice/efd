#pragma once

#include "GridGeometry.hpp"
#include "SolverTraits.hpp"

#include "GhostLayer/IfsfdHandlersBuilder.hpp"
#include "GhostLayer/SfsfdHandlersBuilder.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TScalar,
          int TDimensions,
          int TSolverType,
          int TImmersedBoudnaryType,
          int TDebug>
struct SolverBuilderTraits {};

template <typename TScalar,
          int TDimensions,
          int TImmersedBoudnaryType,
          int TDebug>
struct SolverBuilderTraits<TScalar,
                           TDimensions,
                           0,
                           TImmersedBoudnaryType,
                           TDebug> {
  using SolverTraitsType = SfsfdSolverTraits
                           <UniformGridGeometry<TScalar, TDimensions>,
                            TImmersedBoudnaryType,
                            TDebug,
                            TScalar,
                            TDimensions>;

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = GhostLayer::SfsfdHandlersBuilder
            <SolverTraitsType, TDimension, TDirection>;
};

template <typename TScalar,
          int TDimensions,
          int TImmersedBoudnaryType,
          int TDebug>
struct SolverBuilderTraits<TScalar,
                           TDimensions,
                           1,
                           TImmersedBoudnaryType,
                           TDebug> {
  using SolverTraitsType = IfsfdSolverTraits
                           <UniformGridGeometry<TScalar, TDimensions>,
                            TImmersedBoudnaryType,
                            TDebug,
                            TScalar,
                            TDimensions>;

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = GhostLayer::IfsfdHandlersBuilder
            <SolverTraitsType, TDimension, TDirection>;
};
}
}
