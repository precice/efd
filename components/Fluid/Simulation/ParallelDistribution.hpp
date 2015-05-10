#pragma once

#include <Uni/Helpers/macros>

#include <Eigen/Core>

#include <petscsys.h>

#include <array>

namespace FsiSimulation {
namespace FluidSimulation {
/**
 * TODO:
 * 1. Change name to DomainDecomposition
 */
template <unsigned TDimensions>
class ParallelDistribution {
public:
  enum {
    Dimensions = TDimensions
  };

  using VectorDi =  Eigen::Matrix<int, TDimensions, 1>;

  using Neighbors = std::array<std::array<int, 2>, TDimensions>;

public:
  ParallelDistribution() : mpiCommunicator(PETSC_COMM_WORLD) {}

  ParallelDistribution(ParallelDistribution const& other) = delete;

  ~ParallelDistribution() {}

  ParallelDistribution&
  operator=(ParallelDistribution const& other) = delete;

  int
  rankSize() const;

  void
  initialize(VectorDi const& processorSize_,
             VectorDi const& globalSize_,
             VectorDi const& indent_size);

  int
  getMpiRankFromSubdomainSpatialIndex(VectorDi const& index_) const;

  VectorDi
  getSubdomainSpatialIndexFromMpiRank(int const& rank) const;

  bool
  convertUnindentedSpatialIndexFromGlobalToLocal(
    VectorDi const& global_spatial_index,
    int&            rank,
    VectorDi&       local_spatial_index) const;

  VectorDi
  convertSpatialIndexFromIndentedLocalToUnindentedGlobal(
    VectorDi const& local_spatial_index) const;

  unsigned
  convertUnindentedLocalIndexFromSpatialToSerial(
    VectorDi const& local_spatial_index) const;

  unsigned
  convertUnindentedLocalIndexFromSpatialToSerial(
    int const&      rank,
    VectorDi const& local_spatial_index) const;

  bool
  convertUnindentedGlobalIndexFromSpatialToSerial(
    VectorDi const& global_spatial_index,
    int&            rank,
    unsigned&       serial_index) const;

  VectorDi
  convertSerialIndexToUnindentedLocal(
    int const&            rank,
    unsigned const&       serial_index) const;

  VectorDi
  convertSerialIndexToUnindentedGlobal(
    int const&            rank,
    unsigned const&       serial_index) const;

  std::string
  toString() const;

  VectorDi processorSize;
  VectorDi globalCellSize;

  Uni_PublicProperty(VectorDi, indentSize);
  Uni_PublicProperty(int,      rank);
  VectorDi  localCellSize;
  VectorDi  uniformLocalCellSize;
  VectorDi  lastLocalCellSize;
  VectorDi  index;
  VectorDi  corner;
  Neighbors neighbors;
  MPI_Comm  mpiCommunicator;
};
extern template class ParallelDistribution<2>;
extern template class ParallelDistribution<3>;
}
}
