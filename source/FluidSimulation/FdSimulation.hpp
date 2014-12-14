#ifndef FsiSimulation_FluidSimulation_FdSimulation_hpp
#define FsiSimulation_FluidSimulation_FdSimulation_hpp

#include "CellAccessor.hpp"
#include "Configuration.hpp"
#include "GhostLayer/Handlers.hpp"
#include "Grid.hpp"
#include "LinearSolver.hpp"
#include "ParallelDistribution.hpp"
#include "Parameters.hpp"
#include "Simulation.hpp"
#include "VtkPlot.hpp"
#include "functions.hpp"

#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TGridGeometry,
          typename TMemory,
          typename TScalar,
          int TD>
class FdSimulation : public Simulation {
public:
  typedef Simulation              BaseType;
  typedef typename BaseType::Path Path;
  typedef
    CellAccessor<TGridGeometry, TMemory, TD>
    CellAccessorType;

  typedef typename CellAccessorType::GridGeometryType GridGeometryType;
  typedef typename CellAccessorType::MemoryType       MemoryType;
  typedef typename CellAccessorType::CellType         CellType;
  typedef typename CellType::Velocity                 VelocityType;

  typedef Grid<CellAccessorType, TD>          GridType;
  typedef typename GridType::VectorDi         VectorDiType;
  typedef typename GridGeometryType::VectorDs VectorDsType;
  typedef typename GridType::Factory          CellAccessorFactory;

  typedef Parameters<TScalar, TD> ParametersType;

  typedef ParallelDistribution<TD> ParallelDistributionType;

  typedef
    FluidSimulation::LinearSolver<CellAccessorType, TScalar, TD>
    LinearSolverType;

  typedef typename GhostLayer::Handlers<TD> GhostHandlersType;

  typedef
    VtkPlot<CellAccessorType, TScalar, TD>
    VtkPlotType;

public:
  FdSimulation() {}

  FdSimulation(FdSimulation const& other) = delete;

  ~FdSimulation() {}

  FdSimulation const&
  operator=(FdSimulation const& other) = delete;

  void
  setParameters(VectorDiType const& localSize,
                VectorDsType const& width) {
    _gridGeometry.initialize(
      width,
      _parallelDistribution.globalCellSize,
      _parallelDistribution.corner);

    _memory.allocate(localSize);

    CellAccessorFactory cellAccessorFactory(
      [&] (VectorDiType const& i) {
        return CellAccessorType(i, &_memory, &_gridGeometry);
      });

    _grid.initialize(localSize, cellAccessorFactory);

    logGridInitializationInfo(_grid);
    logParallelTopologyInfo(_parallelDistribution);
  }

  void
  initialize(Path const&        outputDirectory,
             std::string const& fileNamePrefix) {
    _linearSolver.initialize(&_grid,
                             &_parallelDistribution,
                             &_ghostHandler,
                             &_dt);

    for (auto accessor : _grid) {
      accessor.currentCell()->velocity() = VelocityType::Zero();
      accessor.currentCell()->fgh()      = VelocityType::Zero();
      accessor.currentCell()->pressure() = 0.0;
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.velocityInitialization[d][d2]();
      }
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.mpiVelocityExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _maxVelocity    = VelocityType::Zero();
    _dt             = 1;
    _time           = 0;
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

    _dt = TimeStepProcessing<TScalar, TD>::compute(_parameters.re(),
                                                   _parameters.tau(),
                                                   _gridGeometry.minCellWidth(),
                                                   _maxVelocity);
    logInfo("Iteration Number {1}", _iterationCount);
    logInfo("Time step size {1}",   _dt);
    _maxVelocity = VelocityType::Zero();

    for (auto accessor : _grid.innerGrid) {
      typedef FghProcessing<CellAccessorType,
                            ParametersType,
                            TScalar,
                            TD> Fgh;
      Fgh::compute(accessor, _parameters, _dt);
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.fghInitialization[d][d2]();
      }
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.mpiFghExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _linearSolver.solve();

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.mpiPressureExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    for (auto accessor : _grid.innerGrid) {
      typedef VelocityProcessing<CellAccessorType,
                                 TScalar,
                                 TD> VelocityProcessingType;
      VelocityProcessingType::compute(accessor, _dt);
      computeMaxVelocity<CellAccessorType, TScalar, TD>
        (accessor, _maxVelocity);
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.velocityInitialization[d][d2]();
      }
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.mpiVelocityExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _time += _dt;
    ++_iterationCount;

    _plot.plot(_iterationCount, _time, _dt);

    return true;
  }

  MemoryType               _memory;
  GridGeometryType         _gridGeometry;
  GridType                 _grid;
  ParametersType           _parameters;
  ParallelDistributionType _parallelDistribution;
  TScalar                  _dt;
  TScalar                  _time;
  TScalar                  _timeLimit;
  int                      _iterationCount;
  int                      _iterationLimit;
  LinearSolverType         _linearSolver;
  VelocityType             _maxVelocity;
  GhostHandlersType        _ghostHandler;
  VtkPlotType              _plot;
};
}
}
#endif
