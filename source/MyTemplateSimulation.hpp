#ifndef FsiSimulation_MyTemplateSimulation_hpp
#define FsiSimulation_MyTemplateSimulation_hpp

#include "CellAccessor.hpp"
#include "EntryPoint/VtkPlot.hpp"
#include "GhostCellsHandler.hpp"
#include "Grid.hpp"
#include "MySimulation.hpp"
#include "ParallelTopology.hpp"
#include "Parameters.h"
#include "SimulationParameters.hpp"
#include "Solvers/PoissonSolver.hpp"
#include "stencils/mystencils.hpp"

#include <Uni/Logging/macros>

namespace FsiSimulation {
template <typename TGridGeometry,
          typename TMemory,
          typename Scalar,
          int D>
class MyTemplateSimulation : public MySimulation {
public:
  typedef MySimulation        Base;
  typedef typename Base::Path Path;
  typedef
    CellAccessor<TGridGeometry, TMemory, D>
    SpecializedCellAccessor;

  typedef typename SpecializedCellAccessor::GridGeometry GridGeometry;
  typedef typename SpecializedCellAccessor::Memory       Memory;
  typedef typename SpecializedCellAccessor::Cell         Cell;
  typedef typename Cell::Velocity                        Velocity;
  typedef typename Cell::Pressure                        Pressure;

  typedef Grid<SpecializedCellAccessor, D>              SpecializedGrid;
  typedef typename SpecializedGrid::VectorDi            VectorDi;
  typedef typename GridGeometry::VectorDs               VectorDs;
  typedef typename SpecializedGrid::CellAccessorFactory CellAccessorFactory;

  typedef SimulationParameters<Scalar, D> SpecializedSimulationParameters;

  typedef ParallelTopology<D> SpecializedParallelTopology;

  typedef
    Solvers::PoissonSolver<SpecializedCellAccessor, Scalar, D>
    PoissonSolver;

  typedef GhostCellsHandler<D> SpecializedGhostCellsHandler;

  typedef
    EntryPoint::VtkPlot<SpecializedCellAccessor, Scalar, D>
    SpecializedVtkPlot;

public:
  MyTemplateSimulation() {}

  MyTemplateSimulation(MyTemplateSimulation const& other) = delete;

  ~MyTemplateSimulation() {}

  MyTemplateSimulation const&
  operator=(MyTemplateSimulation const& other) = delete;

  void
  setParameters(VectorDi const& localSize,
                VectorDs const& width) {
    _gridGeometry.initialize(
      width,
      _parallelTopology.globalSize,
      _parallelTopology.corner);

    _memory.allocate(localSize);

    CellAccessorFactory cellAccessorFactory(
      [&] (VectorDi const& i) {
        return SpecializedCellAccessor(i, &_memory, &_gridGeometry);
      });

    _grid.initialize(localSize, cellAccessorFactory);

    logGridInitializationInfo(_grid);
    logParallelTopologyInfo(_parallelTopology);
  }

  void
  initialize(Path const&        outputDirectory,
             std::string const& fileNamePrefix) {
    _poissonSolver.initialize(&_grid,
                              &_parallelTopology,
                              &_ghostCellsHandler,
                              &_dt);

    for (auto accessor : _grid) {
      accessor.currentCell()->velocity() = Velocity::Zero();
      accessor.currentCell()->fgh()      = Velocity::Zero();
      accessor.currentCell()->pressure() = 0.0;
    }

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._velocityInitialization[d][d2]();
      }
    }

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._mpiVelocityExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _maxVelocity    = Velocity::Zero();
    _dt             = 1;
    _time           = 0;
    _iterationCount = 0;

    _plot.initialize(&_grid,
                     &_parallelTopology,
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

    _dt = TimeStepProcessing<Scalar, D>::compute(_parameters.re(),
                                                 _parameters.tau(),
                                                 _gridGeometry.minCellWidth(),
                                                 _maxVelocity);
    logInfo("Iteration Number {1}", _iterationCount);
    logInfo("Time step size {1}",   _dt);
    _maxVelocity = Velocity::Zero();

    for (auto accessor : _grid.innerGrid) {
      typedef FghProcessing<SpecializedCellAccessor,
                            SpecializedSimulationParameters,
                            Scalar,
                            D> Fgh;
      Fgh::compute(accessor, _parameters, _dt);
    }

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._fghInitialization[d][d2]();
      }
    }

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._mpiFghExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _poissonSolver.solve();

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._mpiPressureExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    for (auto accessor : _grid.innerGrid) {
      typedef VelocityProcessing<SpecializedCellAccessor, Scalar, D> Velocity;
      Velocity::compute(accessor, _dt);
      computeMaxVelocity<SpecializedCellAccessor, Scalar, D>
        (accessor, _maxVelocity);
    }

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._velocityInitialization[d][d2]();
      }
    }

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler._mpiVelocityExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _time += _dt;
    ++_iterationCount;

    _plot.plot(_iterationCount, _time, _dt);

    return true;
  }

  Memory                          _memory;
  GridGeometry                    _gridGeometry;
  SpecializedGrid                 _grid;
  SpecializedSimulationParameters _parameters;
  SpecializedParallelTopology     _parallelTopology;
  Scalar                          _dt;
  Scalar                          _time;
  Scalar                          _timeLimit;
  int                             _iterationCount;
  int                             _iterationLimit;
  PoissonSolver                   _poissonSolver;
  Velocity                        _maxVelocity;
  SpecializedGhostCellsHandler    _ghostCellsHandler;
  SpecializedVtkPlot              _plot;
};
}
#endif
