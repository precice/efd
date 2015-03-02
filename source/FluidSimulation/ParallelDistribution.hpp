#ifndef FsiSimulation_FluidSimulation_ParallelDistribution_hpp
#define FsiSimulation_FluidSimulation_ParallelDistribution_hpp

#include <Uni/ExecutionControl/exception>
#include <Uni/Logging/format>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <petscsys.h>

#include <array>
#include <cstdlib>

namespace FsiSimulation {
namespace FluidSimulation {
template <int TD>
class ParallelDistribution;

template <int TD>
void
logParallelTopologyInfo(ParallelDistribution<TD> const& topology);

template <int TD>
class ParallelDistribution {
public:
  typedef Eigen::Matrix<int, TD, 1>          VectorDi;
  typedef std::array<std::array<int, 2>, TD> Neighbors;

public:
  ParallelDistribution(): mpiCommunicator(PETSC_COMM_WORLD){}

  ParallelDistribution(ParallelDistribution const& other) = delete;

  ~ParallelDistribution() {}

  ParallelDistribution&
  operator=(ParallelDistribution const& other) = delete;

  void
  initialize(int const&      rank_,
             VectorDi const& processorSize_,
             VectorDi const& globalSize_) {
    int processCount;
    MPI_Comm_size(PETSC_COMM_WORLD, &processCount);

    if (processCount != processorSize_.prod()) {
      throwException("Wrong number of nodes: executed with {1} processors, "
                     "configured with {2} ({3}) processors",
                     processCount,
                     processorSize_.prod(),
                     processorSize_.transpose());
    }

    rank           = rank_;
    processorSize  = processorSize_;
    localCellSize  = globalSize_.cwiseQuotient(processorSize);
    globalCellSize = localCellSize.cwiseProduct(processorSize);

    auto tempDivRank = rank;
    int  tempDivSize = 1;

    for (int i = 0; i < (TD - 1); ++i) {
      tempDivSize = processorSize(i);
      auto tempDiv = std::div(tempDivRank, tempDivSize);
      index(i)    = tempDiv.rem;
      corner(i)   = index(i) * localCellSize(i);
      tempDivRank = tempDiv.quot;
    }
    index(TD - 1)  = tempDivRank;
    corner(TD - 1) = index(TD - 1) * localCellSize(TD - 1);

    for (int d = 0; d < TD; ++d) {
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

    for (int i = 0; i < TD; ++i) {
      if (index_(i) < 0 || index_(i) >= processorSize(i)) {
        return -1;
      }
      result   += index_(i) * tempSize;
      tempSize *= processorSize(i);
    }

    return result;
  }

  VectorDi processorSize;
  VectorDi globalCellSize;

  int       rank;
  VectorDi  localCellSize;
  VectorDi  index;
  VectorDi  corner;
  Neighbors neighbors;
  MPI_Comm  mpiCommunicator;
};

template <int TD>
void
logParallelTopologyInfo(ParallelDistribution<TD> const& topology) {
  std::string neighbors;

  for (int d = 0; d < TD; ++d) {
    neighbors += Uni::Logging::format("{1} ", topology.neighbors[d][0]);
    neighbors += Uni::Logging::format("{1} ", topology.neighbors[d][1]);
  }
  logInfo("ParallelDistribution current rank:     {1}\n"
          "ParallelDistribution processor size:   {2}\n"
          "ParallelDistribution global cell size: {3}\n"
          "ParallelDistribution local cell size:  {4}\n"
          "ParallelDistribution current index:    {5}\n"
          "ParallelDistribution corner:           {6}\n"
          "ParallelDistribution neighbors:        {7}\n",
          topology.rank,
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
