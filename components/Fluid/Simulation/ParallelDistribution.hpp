#ifndef Fluid_Simulation_ParallelDistribution_hpp
#define Fluid_Simulation_ParallelDistribution_hpp

#include <Uni/ExecutionControl/exception>
#include <Uni/Logging/format>
#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <petscsys.h>

#include <array>
#include <cstdlib>

namespace FsiSimulation {
namespace FluidSimulation {
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

  void
  initialize(VectorDi const& processorSize_,
             VectorDi const& globalSize_);

  int
  getRank(VectorDi const& index_) const;

  std::string
  toString() const;

  VectorDi processorSize;
  VectorDi globalCellSize;

  int       rank;
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
#endif
