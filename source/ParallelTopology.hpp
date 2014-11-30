#ifndef FsiSimulation_ParallelTopology_hpp
#define FsiSimulation_ParallelTopology_hpp

#include <Uni/Logging/macros>

#include <Eigen/Core>

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
  typedef Eigen::Matrix<int, D, 1>    VectorDi;
  typedef Eigen::Matrix<int, 2* D, 1> Neighbors;

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
      tempDivSize *= globalSize(i);
      auto tempDiv = std::div(tempDivRank, tempDivSize);
      index(i)    = tempDiv.rem;
      corner(i)   = index(i) * localSize(i);
      tempDivRank = tempDiv.quot;
    }
    index(D - 1)  = tempDivRank;
    corner(D - 1) = index(D - 1) * localSize(D - 1);

    for (int d = 0; d < D; ++d) {
      auto tempIndex = index;
      tempIndex(d)        -= 1;
      neighbors(2 * d)     = getRank(tempIndex);
      tempIndex            = index;
      tempIndex(d)        += 1;
      neighbors(2 * d + 1) = getRank(tempIndex);
    }
  }

  int
  getRank(VectorDi const& index_) {
    auto result   = 0;
    int  tempSize = 1;

    for (int i = 0; i < D; ++i) {
      result   += index_(i) * tempSize;
      tempSize *= globalSize(i);
    }

    if (result < 0 || result >= processorSize.prod()) {
      result = -1;
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
          topology.neighbors.transpose());
}
}

#endif
