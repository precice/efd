#include "ParallelDistribution.hpp"

#include <Uni/ExecutionControl/exception>
#include <Uni/Logging/format>
#include <Uni/Logging/macros>

#include <cstdlib>

namespace FsiSimulation {
namespace FluidSimulation {
template <unsigned T>
int
ParallelDistribution<T>::
rankSize() const {
  return processorSize.prod();
}

template <unsigned T>
void
ParallelDistribution<T>::
initialize(VectorDi const& processorSize_,
           VectorDi const& globalSize_,
           VectorDi const& indent_size) {
  _indentSize = indent_size;

  int processCount;
  MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);
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

  auto tempDivRank = _rank;
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
    neighbors[d][0] = getMpiRankFromSubdomainSpatialIndex(tempIndex);
    tempIndex       = index;
    tempIndex(d)   += 1;
    neighbors[d][1] = getMpiRankFromSubdomainSpatialIndex(tempIndex);
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
getMpiRankFromSubdomainSpatialIndex(VectorDi const& index_) const {
  auto result   = 0;
  int  tempSize = 1;

  for (unsigned d = 0; d < Dimensions; ++d) {
    if (index_(d) < 0 || index_(d) >= processorSize(d)) {
      return -1;
    }
    result   += index_(d) * tempSize;
    tempSize *= processorSize(d);
  }

  return result;
}

template <unsigned T>
typename ParallelDistribution<T>::VectorDi
ParallelDistribution<T>::
getSubdomainSpatialIndexFromMpiRank(int const& rank) const {
  VectorDi spatial_index;

  std::div_t divt;
  divt.quot = rank;

  for (unsigned d = 0; d < Dimensions; ++d) {
    divt             = std::div(divt.quot, processorSize(d));
    spatial_index(d) = divt.rem;
  }

  return spatial_index;
}

template <unsigned T>
bool
ParallelDistribution<T>::
convertUnindentedSpatialIndexFromGlobalToLocal(
  VectorDi const& global_spatial_index,
  int&            rank,
  VectorDi&       local_spatial_index) const {
  VectorDi subdomain_index
    = global_spatial_index.cwiseQuotient(uniformLocalCellSize);

  rank = getMpiRankFromSubdomainSpatialIndex(subdomain_index);

  local_spatial_index
    = global_spatial_index
      - subdomain_index.cwiseProduct(uniformLocalCellSize);

  return rank == _rank;
}

template <unsigned T>
typename ParallelDistribution<T>::VectorDi
ParallelDistribution<T>::
convertSpatialIndexFromIndentedLocalToUnindentedGlobal(
  VectorDi const& local_spatial_index) const {
  return local_spatial_index - _indentSize + corner;
}

template <unsigned T>
unsigned
ParallelDistribution<T>::
convertUnindentedLocalIndexFromSpatialToSerial(
  VectorDi const& local_spatial_index) const {
  auto indented_local_spatial_index = local_spatial_index + _indentSize;

  unsigned local_serial_index = indented_local_spatial_index(Dimensions - 1);

  for (int d = Dimensions - 2; d >= 0; --d) {
    local_serial_index *= localCellSize(d) + 2 * _indentSize(d);
    local_serial_index += indented_local_spatial_index(d);
  }

  return local_serial_index;
}

template <unsigned T>
unsigned
ParallelDistribution<T>::
convertUnindentedLocalIndexFromSpatialToSerial(
  int const&      rank,
  VectorDi const& local_spatial_index) const {
  VectorDi subdomain_spatial_index = getSubdomainSpatialIndexFromMpiRank(rank);

  auto indented_local_spatial_index = local_spatial_index + _indentSize;

  unsigned local_serial_index = indented_local_spatial_index(Dimensions - 1);

  for (int d = Dimensions - 2; d >= 0; --d) {
    if (subdomain_spatial_index(d) == (processorSize(d) - 1)) {
      local_serial_index *= lastLocalCellSize(d) + 2 * _indentSize(d);
      local_serial_index += indented_local_spatial_index(d);
    } else {
      local_serial_index *= uniformLocalCellSize(d) + 2 * _indentSize(d);
      local_serial_index += indented_local_spatial_index(d);
    }
  }

  return local_serial_index;
}

template <unsigned T>
bool
ParallelDistribution<T>::
convertUnindentedGlobalIndexFromSpatialToSerial(
  VectorDi const& global_spatial_index,
  int&            rank,
  unsigned&       serial_index) const {

  assert((global_spatial_index.array() >= VectorDi::Zero().array()).all()
        && (global_spatial_index.array() < globalCellSize.array()).all());
  //

  VectorDi subdomain_index
    = global_spatial_index.cwiseQuotient(uniformLocalCellSize);

  rank = getMpiRankFromSubdomainSpatialIndex(subdomain_index);

  auto local_spatial_index
    = global_spatial_index
      - subdomain_index.cwiseProduct(uniformLocalCellSize);

  bool is_current_domain = rank == _rank;

  if (is_current_domain) {
    serial_index = convertUnindentedLocalIndexFromSpatialToSerial(
      local_spatial_index);
  } else {
    serial_index = convertUnindentedLocalIndexFromSpatialToSerial(
      rank,
      local_spatial_index);
  }

  return is_current_domain;
}

template <unsigned T>
bool
ParallelDistribution<T>::
convertIndentedLocalIndexFromSpatialToSerial(
  VectorDi const& local_spatial_index,
  int&            rank,
  unsigned&       serial_index) const {
  //
  VectorDi global_spatial_index = local_spatial_index - _indentSize + corner;

  return convertUnindentedGlobalIndexFromSpatialToSerial(
         global_spatial_index,
         rank,
         serial_index);
}

template <unsigned T>
typename ParallelDistribution<T>::VectorDi
ParallelDistribution<T>::
convertSerialIndexToUnindentedLocal(
  int const&      rank,
  unsigned const& serial_index) const {
  VectorDi subdomain_spatial_index = getSubdomainSpatialIndexFromMpiRank(rank);

  VectorDi spatial_index;

  std::div_t divt;
  divt.quot = serial_index;

  for (unsigned d = 0; d < Dimensions; ++d) {
    if (subdomain_spatial_index(d) == (processorSize(d) - 1)) {
      divt = std::div(divt.quot, lastLocalCellSize(d) + 2 * _indentSize(d));
    } else {
      divt = std::div(divt.quot, uniformLocalCellSize(d) + 2 * _indentSize(d));
    }
    spatial_index(d) = divt.rem;
  }
  spatial_index -= _indentSize;

  return spatial_index;
}

template <unsigned T>
typename ParallelDistribution<T>::VectorDi
ParallelDistribution<T>::
convertSerialIndexToUnindentedGlobal(
  int const&      rank,
  unsigned const& serial_index) const {
  VectorDi subdomain_spatial_index = getSubdomainSpatialIndexFromMpiRank(rank);

  VectorDi spatial_index;

  std::div_t divt;
  divt.quot = serial_index;

  for (unsigned d = 0; d < Dimensions; ++d) {
    if (subdomain_spatial_index(d) == (processorSize(d) - 1)) {
      divt = std::div(divt.quot, lastLocalCellSize(d) + 2 * _indentSize(d));
    } else {
      divt = std::div(divt.quot, uniformLocalCellSize(d) + 2 * _indentSize(d));
    }
    spatial_index(d) = divt.rem;
  }
  spatial_index -= _indentSize;

  return subdomain_spatial_index.cwiseProduct(uniformLocalCellSize)
         + spatial_index;
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
          this->_rank,
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
