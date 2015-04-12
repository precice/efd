#include "ParallelDistribution.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <unsigned T>
void
ParallelDistribution<T>::
initialize(VectorDi const& processorSize_,
           VectorDi const& globalSize_) {
  int processCount;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &processCount);

  if (processCount != processorSize_.prod()) {
    throwException("Wrong number of nodes: executed with {1} processors, "
                   "configured with {2} ({3}) processors",
                   processCount,
                   processorSize_.prod(),
                   processorSize_.transpose());
  }

  processorSize        = processorSize_;
  uniformLocalCellSize = globalSize_.cwiseQuotient(processorSize);
  localCellSize        = uniformLocalCellSize;
  lastLocalCellSize    = uniformLocalCellSize;
  globalCellSize       = uniformLocalCellSize.cwiseProduct(processorSize);

  auto tempDivRank = rank;
  int  tempDivSize = 1;

  for (unsigned i = 0; i < (Dimensions - 1); ++i) {
    tempDivSize = processorSize(i);
    auto tempDiv = std::div(tempDivRank, tempDivSize);
    index(i)    = tempDiv.rem;
    corner(i)   = index(i) * uniformLocalCellSize(i);
    tempDivRank = tempDiv.quot;
  }
  index(Dimensions - 1)  = tempDivRank;
  corner(Dimensions - 1) = index(Dimensions - 1) * uniformLocalCellSize(
    Dimensions - 1);

  for (unsigned d = 0; d < Dimensions; ++d) {
    auto tempIndex = index;
    tempIndex(d)   -= 1;
    neighbors[d][0] = getRank(tempIndex);
    tempIndex       = index;
    tempIndex(d)   += 1;
    neighbors[d][1] = getRank(tempIndex);
  }

  for (unsigned d = 0; d < Dimensions; ++d) {
    if (globalCellSize(d) != globalSize_(d)) {
      auto diff = globalSize_(d) - globalCellSize(d);

      if (neighbors[d][1] < 0) {
        localCellSize(d) += diff;
      }

      lastLocalCellSize(d) += diff;
      globalCellSize(d)     = globalSize_(d);
    }
  }
}

template <unsigned T>
int
ParallelDistribution<T>::
getRank(VectorDi const& index_) const {
  auto result   = 0;
  int  tempSize = 1;

  for (unsigned i = 0; i < T; ++i) {
    if (index_(i) < 0 || index_(i) >= processorSize(i)) {
      return -1;
    }
    result   += index_(i) * tempSize;
    tempSize *= processorSize(i);
  }

  return result;
}

template <unsigned T>
std::string
ParallelDistribution<T>::
toString() const {
  std::string neighbors;

  for (unsigned d = 0; d < T; ++d) {
    neighbors += Uni::Logging::format("{1} ", this->neighbors[d][0]);
    neighbors += Uni::Logging::format("{1} ", this->neighbors[d][1]);
  }
  logInfo("ParallelDistribution current rank:     {1}\n"
          "ParallelDistribution processor size:   {2}\n"
          "ParallelDistribution global cell size: {3}\n"
          "ParallelDistribution local cell size:  {4}\n"
          "ParallelDistribution current index:    {5}\n"
          "ParallelDistribution corner:           {6}\n"
          "ParallelDistribution neighbors:        {7}\n",
          this->rank,
          this->processorSize.transpose(),
          this->globalCellSize.transpose(),
          this->localCellSize.transpose(),
          this->index.transpose(),
          this->corner.transpose(),
          neighbors);
  return "";
}

template class ParallelDistribution<2>;
template class ParallelDistribution<3>;
}
}
