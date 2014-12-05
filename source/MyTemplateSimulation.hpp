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
  initialize(Path const&        outputDirectory,
             std::string const& fileNamePrefix) {
    VectorDi localSize(_parallelTopology.localSize + 2 * VectorDi::Ones());

    _memory.allocate(localSize);

    CellAccessorFactory cellAccessorFactory(
      [&] (VectorDi const& i) {
        return SpecializedCellAccessor(i, &_memory, &_gridGeometry);
      });

    _grid.initialize(localSize, cellAccessorFactory);

    _poissonSolver.initialize(&_grid,
                              &_parallelTopology,
                              &_ghostCellsHandler,
                              &_dt);

    logGridInitializationInfo(_grid);
    logParallelTopologyInfo(_parallelTopology);

    for (auto accessor : _grid) {
      accessor.currentCell()->velocity() = Velocity::Zero();
      accessor.currentCell()->fgh()      = Velocity::Zero();
      accessor.currentCell()->pressure() = 0.0;
    }

    //for (int d = 0; d < D; ++d) {
    //  for (int d2 = 0; d2 < 2; ++d2) {
    //    _ghostCellsHandler._velocityInitialization[d][d2]();
    //  }
    //}

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
    // logInfo("Maximumal Velocity {1}", _maxVelocity.transpose());
    // logInfo("re {1}",                 _parameters.re());
    // logInfo("gamma {1}",              _parameters.gamma());
    // logInfo("tau {1}",                _parameters.tau());
    // logInfo("g {1}",                  _parameters.g().transpose());
    // logInfo("min cell width {1}",
    // _gridGeometry.minCellWidth().transpose());
    //
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

    _poissonSolver.solve();

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
