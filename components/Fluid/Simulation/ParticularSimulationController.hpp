#pragma once

#include "SimulationController.hpp"
#include "SolverTraits.hpp"

#include <Uni/Firewall/Implementation>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class ParticularSimulationControllerImplementation;

class IterationResultWriter;
class Reporter;

template <typename TSolverTraits>
class ParticularSimulationController : public SimulationController {
private:
  using Implementation
          = ParticularSimulationControllerImplementation<TSolverTraits>;

public:
  using BaseType = SimulationController;

  using Path = typename BaseType::Path;

  using SolverTraitsType = TSolverTraits;

  using SolverType = typename SolverTraitsType::SolverType;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridGeometryType = typename SolverTraitsType::GridGeometryType;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

public:
  ParticularSimulationController();

  ParticularSimulationController(
    ParticularSimulationController const& other) = delete;

  ~ParticularSimulationController();

  ParticularSimulationController const&
  operator=(ParticularSimulationController const& other) = delete;

  void
  setIterationResultWriter(IterationResultWriter* iteration_result_writer);

  SolverType const*
  solver() const;

  SolverType*
  solver();

  long double const&
  time() const;

  long double const&
  timeLimit() const;

  long double&
  timeLimit();

  long double const&
  lastPlotTimeStamp() const;

  long double&
  lastPlotTimeStamp();

  ScalarType const&
  plotInterval() const;

  ScalarType&
  plotInterval();

  unsigned long long const&
  iterationNumber() const;

  unsigned long long const&
  iterationLimit() const;

  unsigned long long&
  iterationLimit();

  void
  initialize(precice::SolverInterface* preciceInteface,
             Reporter*                 reporter,
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix);

  bool
  iterate();

private:
  Uni_Firewall_IMPLEMENTATION_LINK(
    ParticularSimulationControllerImplementation<TSolverTraits> );
};

Fluid_DeclareExternTemplates(ParticularSimulationController);
}
}
