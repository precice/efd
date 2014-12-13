#ifndef FsiSimulation_MyTemplateSimulation_hpp
#define FsiSimulation_MyTemplateSimulation_hpp

#include "ParallelDistribution.hpp"
#include "CellAccessor.hpp"
#include "VtkPlot.hpp"
#include "GhostLayer/Handlers.hpp"
#include "Grid.hpp"
#include "Simulation.hpp"
#include "Configuration.hpp"
#include "SimulationParameters.hpp"
#include "LinearSolver.hpp"
#include "functions.hpp"

#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template<typename TGridGeometry,
typename TMemory,
typename Scalar,
int TD>
class FdSimulation : public Simulation {
public:
  typedef Simulation BaseType;
  typedef typename BaseType::Path Path;
  typedef
  CellAccessor<TGridGeometry, TMemory, TD>
  CellAccessorType;

  typedef typename CellAccessorType::GridGeometryType GridGeometryType;
  typedef typename CellAccessorType::MemoryType MemoryType;
  typedef typename CellAccessorType::CellType CellType;
  typedef typename CellType::Velocity VelocityType;

  typedef Grid<CellAccessorType, TD> GridType;
  typedef typename GridType::VectorDi VectorDi;
  typedef typename GridGeometryType::VectorDs VectorDs;
  typedef typename GridType::CellAccessorFactory CellAccessorFactory;

  typedef SimulationParameters<Scalar, TD> SpecializedSimulationParameters;

  typedef ParallelDistribution<TD> SpecializedParallelTopology;

  typedef
  FluidSimulation::LinearSolver<CellAccessorType, Scalar, TD>
  PoissonSolver;

  typedef typename GhostLayer::Handlers<TD> SpecializedGhostCellsHandler;

  typedef
  VtkPlot<CellAccessorType, Scalar, TD>
  SpecializedVtkPlot;

public:
  FdSimulation() {
  }

  FdSimulation(FdSimulation const& other) = delete;

  ~FdSimulation() {
  }

  FdSimulation const&
  operator=(FdSimulation const& other) = delete;

  void
  setParameters(VectorDi const& localSize,
  VectorDs const& width) {
    _gridGeometry.initialize(
    width,
    _parallelDistribution.globalCellSize,
    _parallelDistribution.corner);

    _memory.allocate(localSize);

    CellAccessorFactory cellAccessorFactory(
    [&](VectorDi const& i) {
      return CellAccessorType(i, &_memory, &_gridGeometry);
    });

    _grid.initialize(localSize, cellAccessorFactory);

    logGridInitializationInfo(_grid);
    logParallelTopologyInfo(_parallelDistribution);
  }

  void
  initialize(Path const& outputDirectory,
  std::string const& fileNamePrefix) {
    _poissonSolver.initialize(&_grid,
    &_parallelDistribution,
    &_ghostCellsHandler,
    &_dt);

    for (auto accessor : _grid) {
      accessor.currentCell()->velocity() = VelocityType::Zero();
      accessor.currentCell()->fgh() = VelocityType::Zero();
      accessor.currentCell()->pressure() = 0.0;
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._velocityInitialization[d][d2]();
      }
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._mpiVelocityExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _maxVelocity = VelocityType::Zero();
    _dt = 1;
    _time = 0;
    _iterationCount = 0;

    _plot.initialize(&_grid,
    &_parallelDistribution,
    &_gridGeometry,
    outputDirectory,
    fileNamePrefix);

    _plot.plot(_iterationCount, _time, _dt);
  }

  bool
  iterate() {
    if (_time >= _timeLimit ||
    _iterationCount >= _iterationLimit) {
      return false;
    }

    _dt = TimeStepProcessing<Scalar, TD>::compute(_parameters.re(),
    _parameters.tau(),
    _gridGeometry.minCellWidth(),
    _maxVelocity);
    logInfo("Iteration Number {1}", _iterationCount);
    logInfo("Time step size {1}", _dt);
    _maxVelocity = VelocityType::Zero();

    for (auto accessor : _grid.innerGrid) {
      typedef FghProcessing<CellAccessorType,
      SpecializedSimulationParameters,
      Scalar,
      TD> Fgh;
      Fgh::compute(accessor, _parameters, _dt);
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._fghInitialization[d][d2]();
      }
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._mpiFghExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _poissonSolver.solve();

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._mpiPressureExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    for (auto accessor : _grid.innerGrid) {
      typedef VelocityProcessing<CellAccessorType, Scalar, TD> VelocityProcessingType;
      VelocityProcessingType::compute(accessor, _dt);
      computeMaxVelocity<CellAccessorType, Scalar, TD>
      (accessor, _maxVelocity);
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._velocityInitialization[d][d2]();
      }
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._mpiVelocityExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _time += _dt;
    ++_iterationCount;

    _plot.plot(_iterationCount, _time, _dt);

    return true;
  }

  MemoryType _memory;
  GridGeometryType _gridGeometry;
  GridType _grid;
  SpecializedSimulationParameters _parameters;
  SpecializedParallelTopology _parallelDistribution;
  Scalar _dt;
  Scalar _time;
  Scalar _timeLimit;
  int _iterationCount;
  int _iterationLimit;
  PoissonSolver _poissonSolver;
  VelocityType _maxVelocity;
  SpecializedGhostCellsHandler _ghostCellsHandler;
  SpecializedVtkPlot _plot;
};
}
}
#endif
