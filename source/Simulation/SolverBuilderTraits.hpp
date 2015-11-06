#pragma once

#include "SolverTraits.hpp"

namespace Fluid {
namespace Simulation {
namespace GhostLayer {
template <typename T, int U, int D>
class SfsfdHandlersBuilder;
template <typename T, int U, int D>
class IfsfdHandlersBuilder;
}
template <int TSolverId,
          int TDebug,
          typename TScalar,
          int TDimensions>
struct SolverBuilderTraits {};

template <int TDebug,
          typename TScalar,
          int TDimensions>
struct SolverBuilderTraits<0,
                           TDebug,
                           TScalar,
                           TDimensions> {
  using SolverTraitsType = SfsfdSolverTraits
                           < TDebug,
                            TScalar,
                            TDimensions>;

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = GhostLayer::SfsfdHandlersBuilder
            <SolverTraitsType, TDimension, TDirection>;
};

template <int TDebug,
          typename TScalar,
          int TDimensions>
struct SolverBuilderTraits<1,
                           TDebug,
                           TScalar,
                           TDimensions> {
  using SolverTraitsType = IfsfdSolverTraits
                           <TDebug,
                            TScalar,
                            TDimensions>;
  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = GhostLayer::IfsfdHandlersBuilder
            <SolverTraitsType, TDimension, TDirection>;
};
}
}
