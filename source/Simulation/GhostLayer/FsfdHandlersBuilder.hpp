#pragma once

#include "Simulation/Grid.hpp"

#include "Simulation/SolverTraits.hpp"

namespace Fluid {
namespace Simulation {
class Configuration;
namespace GhostLayer {
namespace MpiExchange {
template <typename TDataElement,
          unsigned TDataElementSize,
          typename TGrid,
          TDataElement* (*)(typename TGrid::CellAccessor const&),
          void (*)(typename TGrid::CellAccessor const&,
                   int const&,
                   TDataElement const&),
          int TDimension,
          int TDirection>
class Handler;
}

template <typename THandlersBuilderTraits>
class FsfdHandlersBuilder {
public:
  using SolverTraitsType = typename THandlersBuilderTraits::SolverTraitsType;

  using SolverType = typename SolverTraitsType::SolverType;

  enum {
    Dimension  = THandlersBuilderTraits::Dimension,
    Direction  = THandlersBuilderTraits::Direction,
    Dimensions = SolverTraitsType::Dimensions
  };

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  template <int TDimension>
  inline static ScalarType*
  getFgh(CellAccessorType const& accessor) {
    return &accessor.fgh(TDimension);
  }

  template <int TDimension>
  inline static void
  setFgh(CellAccessorType const& accessor,
         int const&              index,
         ScalarType const&       value) {
    ((void)index);
    accessor.fgh(TDimension) = value;
  }

  inline static ScalarType*
  getVelocity(CellAccessorType const& accessor) {
    return accessor.velocity().data();
  }

  inline static void
  setVelocity(CellAccessorType const& accessor,
              int const&              index,
              ScalarType const&       value) {
    accessor.velocity(index) = value;
  }

  inline static int*
  getLocations(CellAccessorType const& accessor) {
    return accessor.positionInRespectToGeometry().data();
  }

  inline static void
  setLocations(CellAccessorType const& accessor,
               int const&              index,
               int const&              value) {
    accessor.positionInRespectToGeometry(index) = value;
  }

public:
  FsfdHandlersBuilder(Configuration*, SolverType*);

  void
  setAsInput();

  void
  setAsParabolicInput();

  void
  setAsOutput();

  void
  setAsMpiExchange();

protected:
  Configuration* _configuration;
  SolverType*    _solver;
};
}
}
}
