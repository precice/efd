#pragma once

#include "SolverTraits.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <typename T, int U, int D>
class SfsfdHandlersBuilder;
template <typename T, int U, int D>
class IfsfdHandlersBuilder;
}
template <int TSolverId,
          int TImmersedBoudnaryType,
          int TDebug,
          typename TScalar,
          int TDimensions>
struct SolverBuilderTraits {};

template <int TImmersedBoudnaryType,
          int TDebug,
          typename TScalar,
          int TDimensions>
struct SolverBuilderTraits<0,
                           TImmersedBoudnaryType,
                           TDebug,
                           TScalar,
                           TDimensions> {
  using SolverTraitsType = SfsfdSolverTraits
                           <TImmersedBoudnaryType,
                            TDebug,
                            TScalar,
                            TDimensions>;

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = GhostLayer::SfsfdHandlersBuilder
            <SolverTraitsType, TDimension, TDirection>;
};

template <int TImmersedBoudnaryType,
          int TDebug,
          typename TScalar,
          int TDimensions>
struct SolverBuilderTraits<1,
                           TImmersedBoudnaryType,
                           TDebug,
                           TScalar,
                           TDimensions> {
  using SolverTraitsType = IfsfdSolverTraits
                           <TImmersedBoudnaryType,
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
