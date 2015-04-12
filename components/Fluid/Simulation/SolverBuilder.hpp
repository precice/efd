#pragma once

#include "SolverBuilderTraits.hpp"

#include "Configuration.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
class Configuration;
template <typename TSolverBuilderTraits>
class SolverBuilder {
private:
  using SolverBuilderTraitsType = TSolverBuilderTraits;

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = typename SolverBuilderTraitsType::template GhostHandlersBuilderType
            <TDimension, TDirection>;

  using SolverTraitsType = typename SolverBuilderTraitsType::SolverTraitsType;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

public:
  using SolverType
          = typename SolverTraitsType::SolverType;

public:
  SolverBuilder(Configuration*, SolverType*);

  void
  setLeftWallAs(WallEnum type);

  void
  setRightWallAs(WallEnum type);

  void
  setBottomWallAs(WallEnum type);

  void
  setTopWallAs(WallEnum type);

  void
  setBackWallAs(WallEnum type);

  void
  setFrontWallAs(WallEnum type);

private:
  Configuration* _configuration;
  SolverType*    _solver;
};

extern template class SolverBuilder
  < SolverBuilderTraits < 0, 0, 0, double, 2 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 0, 0, 1, double, 2 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 0, 1, 0, double, 2 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 0, 1, 1, double, 2 >>;

extern template class SolverBuilder
  < SolverBuilderTraits < 0, 0, 0, double, 3 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 0, 0, 1, double, 3 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 0, 1, 0, double, 3 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 0, 1, 1, double, 3 >>;

extern template class SolverBuilder
  < SolverBuilderTraits < 1, 0, 0, double, 2 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 1, 0, 1, double, 2 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 1, 1, 0, double, 2 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 1, 1, 1, double, 2 >>;

extern template class SolverBuilder
  < SolverBuilderTraits < 1, 0, 0, double, 3 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 1, 0, 1, double, 3 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 1, 1, 0, double, 3 >>;
extern template class SolverBuilder
  < SolverBuilderTraits < 1, 1, 1, double, 3 >>;
}
}
