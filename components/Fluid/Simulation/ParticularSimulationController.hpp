#pragma once

#include "FsfdSolver.hpp"
#include "IterationResultWriter.hpp"
#include "Reporter.hpp"
#include "SimulationController.hpp"

#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class ParticularSimulationController : public SimulationController {
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
  ParticularSimulationController() {}

  ParticularSimulationController(
    ParticularSimulationController const& other) = delete;

  ~ParticularSimulationController() {}

  ParticularSimulationController const&
  operator=(ParticularSimulationController const& other) = delete;

  void
  setIterationResultWriter(IterationResultWriter* iteration_result_writer) {
    _resultWriter.reset(iteration_result_writer);
  }

  SolverType const*
  solver() const {
    return &_solver;
  }

  SolverType*
  solver() {
    return &_solver;
  }

  long double const&
  time() const {
    return _solver.memory()->time();
  }

  long double const&
  timeLimit() const {
    return _timeLimit;
  }

  long double&
  timeLimit() {
    return _timeLimit;
  }

  long double const&
  lastPlotTimeStamp() const {
    return _lastPlotTimeStamp;
  }

  long double&
  lastPlotTimeStamp() {
    return _lastPlotTimeStamp;
  }

  ScalarType const&
  plotInterval() const {
    return _plotInterval;
  }

  ScalarType&
  plotInterval() {
    return _plotInterval;
  }

  unsigned long long const&
  iterationNumber() const {
    return _solver.memory()->iterationNumber();
  }

  unsigned long long const&
  iterationLimit() const {
    return _iterationLimit;
  }

  unsigned long long&
  iterationLimit() {
    return _iterationLimit;
  }

  void
  initialize(precice::SolverInterface* preciceInteface,
             Reporter*                 reporter,
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix) {
    // logInfo("Enter Initialize");
    _reporter = reporter;

    _reporter->setAt("Name", fileNamePrefix);

    // logInfo("Solver is being initialized");
    _solver.initialize(preciceInteface, reporter);

    // logInfo("Result Writer is being set its destination");
    _resultWriter->setDestination(outputDirectory, fileNamePrefix);

    _lastPlotTimeStamp = 0.0;
    _lastTime          = -1.0;

    if (_plotInterval >= 0) {
      _resultWriter->writeGeometry();
      _resultWriter->writeAttributes();
    }
  }

  bool
  iterate() {
    if (std::abs(_lastTime - _solver.memory()->time())
        <= std::numeric_limits<long double>::epsilon()) {
      return false;
    }

    if (_timeLimit > 0
        && _solver.memory()->time() >= _timeLimit) {
      return false;
    }

    if (_iterationLimit > 0
        && _solver.memory()->iterationNumber() >= _iterationLimit) {
      return false;
    }
    _lastTime = _solver.memory()->time();

    _solver.iterate();

    if (_plotInterval >= 0) {
      if ((_solver.memory()->time() - _lastPlotTimeStamp) > _plotInterval) {
        _lastPlotTimeStamp = _solver.memory()->time();

        _resultWriter->writeAttributes();
      }
    }

    return true;
  }

private:
  SolverType                             _solver;
  std::unique_ptr<IterationResultWriter> _resultWriter;

  long double        _lastTime;
  long double        _timeLimit;
  long double        _lastPlotTimeStamp;
  ScalarType         _plotInterval;
  unsigned long long _iterationLimit;

  Reporter* _reporter;
};
}
}
