#pragma once

#include "GridGeometry.hpp"
#include "SfsfdGhostHandlersBuilder.hpp"
#include "SfsfdSolver.hpp"
#include "VtkOutput/VtkWriter.hpp"
#include "SolverTraits.hpp"
#include "XdmfHdf5Output/XdmfHdf5Writer.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TScalar,
          int TDimensions,
          int TSolverType,
          int TOutputWriter>
struct SolverBuilderTraits {};

template <typename TScalar, int TDimensions>
struct SolverBuilderTraits<TScalar, TDimensions, 0, 0> {


  using GridGeometryType = UniformGridGeometry<TScalar, TDimensions>;

  using SolverTraitsType = SolverTraits<GridGeometryType,
                                        TScalar,
                                        TDimensions>;


  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  using GhostHandlersType = typename GhostLayer::SfsfdHandlers<Dimensions>;

  using OutputWriterType = XdmfHdf5Output::XdmfHdf5Writer<MemoryType>;

  using SolverType = SfsfdSolver<GridGeometryType,
                                 OutputWriterType,
                                 ScalarType,
                                 Dimensions>;

  template <int TDimension, int TDirection>
  using GhostHandlersBuilderType
          = SfsfdGhostHandlersBuilder<SolverType,
                                      TDimension,
                                      TDirection>;
};
}
}
