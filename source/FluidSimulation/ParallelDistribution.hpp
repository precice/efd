#ifndef FsiSimulation_ParallelTopology_hpp
#define FsiSimulation_ParallelTopology_hpp

#include <Uni/Logging/format>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <array>
#include <cstdlib>

namespace FsiSimulation {
namespace FluidSimulation {
template<int D>
class ParallelDistribution;

template<int D>
void
logParallelTopologyInfo(ParallelDistribution<D> const& topology);

template<int D>
class ParallelDistribution {
public:
  typedef Eigen::Matrix<int, D, 1> VectorDi;
  typedef std::array<std::array<int, 2>, D> Neighbors;

public:
  ParallelDistribution() {
  }

  ParallelDistribution(ParallelDistribution const& other) = delete;

  ~ParallelDistribution() {
  }

  ParallelDistribution&
  operator=(ParallelDistribution const& other) = delete;

  void
  initialize(int const& rank_,
  VectorDi const& processorSize_,
  VectorDi const& globalSize_) {
    currentRank = rank_;
    processorSize = processorSize_;
    localCellSize = globalSize_.cwiseQuotient(processorSize);
    globalCellSize = localCellSize.cwiseProduct(processorSize);

    auto tempDivRank = currentRank;
    int tempDivSize = 1;

    for (int i = 0; i < (D - 1); ++i) {
      tempDivSize = processorSize(i);
      auto tempDiv = std::div(tempDivRank, tempDivSize);
      index(i) = tempDiv.rem;
      corner(i) = index(i) * localCellSize(i);
      tempDivRank = tempDiv.quot;
    }
    index(D - 1) = tempDivRank;
    corner(D - 1) = index(D - 1) * localCellSize(D - 1);

    for (int d = 0; d < D; ++d) {
      auto tempIndex = index;
      tempIndex(d) -= 1;
      neighbors[d][0] = getRank(tempIndex);
      tempIndex = index;
      tempIndex(d) += 1;
      neighbors[d][1] = getRank(tempIndex);
    }
  }

  int
  getRank(VectorDi const& index_) {
    auto result = 0;
    int tempSize = 1;

    for (int i = 0; i < D; ++i) {
      if (index_(i) < 0 || index_(i) >= processorSize(i)) {
        return -1;
      }
      result += index_(i) * tempSize;
      tempSize *= processorSize(i);
    }

    return result;
  }

  VectorDi processorSize;
  VectorDi globalCellSize;

  int currentRank;
  VectorDi localCellSize;
  VectorDi index;
  VectorDi corner;
  Neighbors neighbors;
};

template<int D>
void
logParallelTopologyInfo(ParallelDistribution<D> const& topology) {
  std::string neighbors;

  for (int d = 0; d < D; ++d) {
    neighbors += Uni::Logging::format("{1} ", topology.neighbors[d][0]);
    neighbors += Uni::Logging::format("{1} ", topology.neighbors[d][1]);
  }
  logInfo("ParallelDistribution current rank:     {1}\n"
  "ParallelDistribution processor size:   {2}\n"
  "ParallelDistribution global cell size: {3}\n"
  "ParallelDistribution local cell size:  {4}\n"
  "ParallelDistribution index:            {5}\n"
  "ParallelDistribution corner:           {6}\n"
  "ParallelDistribution neighbors:        {7}\n",
  topology.currentRank,
  topology.processorSize.transpose(),
  topology.globalCellSize.transpose(),
  topology.localCellSize.transpose(),
  topology.index.transpose(),
  topology.corner.transpose(),
  neighbors);
}
}
}

#endif
