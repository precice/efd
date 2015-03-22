#pragma once

#include "SimulationController.hpp"

#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits,
          typename ResultWriter>
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

  using ScalarType = typename SolverType::ScalarType;

public:
  ParticularSimulationController() {}

  ParticularSimulationController(
    ParticularSimulationController const& other) = delete;

  ~ParticularSimulationController() {}

  ParticularSimulationController const&
  operator=(ParticularSimulationController const& other) = delete;

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
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix) {
    _solver.initialize(preciceInteface);

    _resultWriter.initialize(_solver.memory(), outputDirectory, fileNamePrefix);

    _lastPlotTimeStamp = 0.0;

    if (_plotInterval >= 0) {
      _resultWriter.writeGeometry();
      _resultWriter.writeAttributes();
    }
  }

  bool
  iterate() {
    if (_timeLimit > 0
        && _solver.memory()->time() >= _timeLimit) {
      return false;
    }

    if (_iterationLimit > 0
        && _solver.memory()->iterationNumber() >= _iterationLimit) {
      return false;
    }

    _solver.iterate();

    if (_plotInterval >= 0) {
      if ((_solver.memory()->time() - _lastPlotTimeStamp) > _plotInterval) {
        _lastPlotTimeStamp = _solver.memory()->time();

        _resultWriter.writeAttributes();
      }
    }

    return true;
  }

private:
  SolverType _solver;

  long double        _timeLimit;
  long double        _lastPlotTimeStamp;
  ScalarType         _plotInterval;
  unsigned long long _iterationLimit;

  ResultWriter _resultWriter;
};
}
}
