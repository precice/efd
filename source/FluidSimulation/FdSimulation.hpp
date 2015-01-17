#ifndef FsiSimulation_FluidSimulation_FdSimulation_hpp
#define FsiSimulation_FluidSimulation_FdSimulation_hpp

#include "Configuration.hpp"
#include "GhostLayer/Handlers.hpp"
#include "Grid.hpp"
#include "ImmersedBoundary/BodyForce/functions.hpp"
#include "ImmersedBoundary/Fadlun/functions.hpp"
#include "ImmersedBoundary/FeedbackForcing/functions.hpp"
#include "LinearSolver.hpp"
#include "ParallelDistribution.hpp"
#include "Parameters.hpp"
#include "Private/mpigenerics.hpp"
#include "Simulation.hpp"
#include "VtkPlot.hpp"
#include "functions.hpp"

#include <Uni/Logging/macros>
#include <Uni/Stopwatch>

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

  typedef Grid<TMemory, TGridGeometry, TD>            GridType;
  typedef typename GridType::CellAccessorType         CellAccessorType;
  typedef typename CellAccessorType::MemoryType       MemoryType;
  typedef typename CellAccessorType::CellType         CellType;
  typedef typename CellAccessorType::GridGeometryType GridGeometryType;
  typedef typename CellType::Velocity                 VelocityType;
  typedef typename GridType::VectorDi                 VectorDiType;
  typedef typename GridGeometryType::VectorDs         VectorDsType;
  typedef typename GridType::Factory                  CellAccessorFactory;

  typedef Parameters<TScalar, TD> ParametersType;

  typedef ParallelDistribution<TD> ParallelDistributionType;

  typedef
    FluidSimulation::LinearSolver<TMemory, TGridGeometry, TScalar, TD>
    LinearSolverType;

  typedef typename GhostLayer::Handlers<TD> GhostHandlersType;

  typedef
    VtkPlot<TMemory, TGridGeometry, TScalar, TD>
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
        return CellAccessorType(i, &_memory, &_grid.innerGrid, &_gridGeometry);
      });

    _grid.initialize(localSize, cellAccessorFactory);

    // logGridInitializationInfo(_grid);
    logParallelTopologyInfo(_parallelDistribution);
  }

  void
  initialize(precice::SolverInterface* preciceInteface,
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix) {
    _preciceInteface = preciceInteface;
    _linearSolver.initialize(&_grid,
                             &_parallelDistribution,
                             &_ghostHandler,
                             &_dt);

    for (auto const& accessor : _grid) {
      accessor.currentCell()->velocity() = VelocityType::Zero();
      accessor.currentCell()->fgh()      = VelocityType::Zero();
      accessor.currentCell()->pressure() = 0.0;
      ImmersedBoundary::Fadlun::template computeDistances<CellAccessorType, TD>(
        accessor,
        _preciceInteface);
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

    _maxVelocity       = VelocityType::Zero();
    _dt                = 1.0;
    _time              = 0.0;
    _lastPlotTimeStamp = 0.0;
    _iterationCount    = 0;

    if (_plotInterval >= 0) {
      _plot.initialize(&_grid,
                       &_parallelDistribution,
                       &_gridGeometry,
                       outputDirectory,
                       fileNamePrefix);

      _plot.plot(_iterationCount, _time, _dt);
    }
  }

  bool
  iterate() {
    if (_timeLimit > 0 && _time >= _timeLimit) {
      return false;
    }

    if (_iterationLimit > 0 && _iterationCount >= _iterationLimit) {
      return false;
    }

    _dt = TimeStepProcessing<TScalar, TD>::compute(_parameters.re(),
                                                   _parameters.tau(),
                                                   _gridGeometry.minCellWidth(),
                                                   _maxVelocity);
    _maxVelocity = VelocityType::Zero();

    if (_parallelDistribution.rank == 0) {
      logInfo("N = {1}; t = {2}", _iterationCount, _time);
    }
    // logInfo("Time step size {1}",   _dt);

    // for (auto accessor : _grid.innerGrid) {
    // for (int d = 0; d < TD; ++d) {
    // if (!ImmersedBoundary::Fadlun::treatBoundary
    // <CellAccessorType, TScalar, TD>(
    // accessor,
    // _preciceInteface,
    // _dt,
    // d)) {
    ////
    // } else {
    // computeMaxVelocity<CellAccessorType, TScalar, TD>
    // (accessor, _maxVelocity);
    // }
    // }
    // }

    for (auto accessor : _grid.innerGrid) {
      typedef FghProcessing<CellAccessorType,
                            ParametersType,
                            TScalar,
                            TD> Fgh;
      Fgh::compute(accessor, _parameters, _dt);

      for (int d = 0; d < TD; ++d) {
        if (ImmersedBoundary::FeedbackForcing
            ::template treatBoundary<CellAccessorType, TD>(accessor,
                                                           _parameters.alpha(),
                                                           d)) {
          //
        }
      }
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

    for (auto const& accessor : _grid.innerGrid) {
      typedef VelocityProcessing<CellAccessorType,
                                 TScalar,
                                 TD> VelocityProcessingType;

      for (int d = 0; d < TD; ++d) {
        VelocityProcessingType::compute(accessor, d, _dt);
      }
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

    VectorDsType force = VectorDsType::Zero();

    for (auto const& accessor : _grid.innerGrid) {
      ImmersedBoundary::BodyForce::
      template computeCellForce<CellAccessorType, TD>(
        accessor,
        _parameters.re(),
        force);
    }

    Private::mpiAllReduce<TScalar>(MPI_IN_PLACE,
                                   force.data(),
                                   TD,
                                   MPI_SUM,
                                   PETSC_COMM_WORLD);

    _time += _dt;
    ++_iterationCount;

    if (_plotInterval >= 0) {
      if ((_time - _lastPlotTimeStamp) > _plotInterval) {
        _lastPlotTimeStamp = _time;
        _plot.plot(_iterationCount, _time, _dt);
      }
    }

    force *= 2 / (0.3 * 0.3 * 0.1);

    if (_parallelDistribution.rank == 0) {
      logInfo("{1}", force.transpose());
    }

    return true;
  }

  MemoryType                _memory;
  GridGeometryType          _gridGeometry;
  GridType                  _grid;
  ParametersType            _parameters;
  ParallelDistributionType  _parallelDistribution;
  precice::SolverInterface* _preciceInteface;
  TScalar                   _dt;
  TScalar                   _time;
  TScalar                   _timeLimit;
  TScalar                   _lastPlotTimeStamp;
  TScalar                   _plotInterval;
  int                       _iterationCount;
  int                       _iterationLimit;
  LinearSolverType          _linearSolver;
  VelocityType              _maxVelocity;
  GhostHandlersType         _ghostHandler;
  VtkPlotType               _plot;
};
}
}
#endif
