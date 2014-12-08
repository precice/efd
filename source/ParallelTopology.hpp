#ifndef FsiSimulation_ParallelTopology_hpp
#define FsiSimulation_ParallelTopology_hpp

#include <Uni/Logging/format>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <array>
#include <cstdlib>

namespace FsiSimulation {
template <int D>
class ParallelTopology;
template <int D>
void
logParallelTopologyInfo(ParallelTopology<D> const& topology);
template <int D>
class ParallelTopology {
public:
  typedef Eigen::Matrix<int, D, 1>          VectorDi;
  typedef std::array<std::array<int, 2>, D> Neighbors;

public:
  ParallelTopology() {}

  ParallelTopology(ParallelTopology const& other) = delete;

  ~ParallelTopology() {}

  ParallelTopology&
  operator=(ParallelTopology const& other) = delete;

  void
  initialize(int const&      rank_,
             VectorDi const& processorSize_,
             VectorDi const& globalSize_) {
    rank          = rank_;
    processorSize = processorSize_;
    localSize     = globalSize_.cwiseQuotient(processorSize);
    globalSize    = localSize.cwiseProduct(processorSize);

    auto tempDivRank = rank;
    int  tempDivSize = 1;

    for (int i = 0; i < (D - 1); ++i) {
      tempDivSize = processorSize(i);
      auto tempDiv = std::div(tempDivRank, tempDivSize);
      index(i)    = tempDiv.rem;
      corner(i)   = index(i) * localSize(i);
      tempDivRank = tempDiv.quot;
    }
    index(D - 1)  = tempDivRank;
    corner(D - 1) = index(D - 1) * localSize(D - 1);

    for (int d = 0; d < D; ++d) {
      auto tempIndex = index;
      tempIndex(d)   -= 1;
      neighbors[d][0] = getRank(tempIndex);
      tempIndex       = index;
      tempIndex(d)   += 1;
      neighbors[d][1] = getRank(tempIndex);
    }
  }

  int
  getRank(VectorDi const& index_) {
    auto result   = 0;
    int  tempSize = 1;

    for (int i = 0; i < D; ++i) {
      if (index_(i) < 0 || index_(i) >= processorSize(i)) {
        return -1;
      }
      result   += index_(i) * tempSize;
      tempSize *= processorSize(i);
    }

    return result;
  }

  VectorDi processorSize;
  VectorDi globalSize;

  int       rank;
  VectorDi  localSize;
  VectorDi  index;
  VectorDi  corner;
  Neighbors neighbors;
};

template <int D>
void
logParallelTopologyInfo(ParallelTopology<D> const& topology) {
  std::string neighbors;

  for (int d = 0; d < D; ++d) {
    neighbors += Uni::Logging::format("{1} ", topology.neighbors[d][0]);
    neighbors += Uni::Logging::format("{1} ", topology.neighbors[d][1]);
  }
  logInfo("Topology rank: {1}\n"
          "Topology processor size: {2}\n"
          "Topology global size: {3}\n"
          "Topology local size: {4}\n"
          "Topology index: {5}\n"
          "Topology corner: {6}\n"
          "Topology neighbors: {7}\n",
          topology.rank,
          topology.processorSize.transpose(),
          topology.globalSize.transpose(),
          topology.localSize.transpose(),
          topology.index.transpose(),
          topology.corner.transpose(),
          neighbors);
}
}

#endif
