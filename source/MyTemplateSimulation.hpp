#ifndef FsiSimulation_MyTemplateSimulation_hpp
#define FsiSimulation_MyTemplateSimulation_hpp

#include "CellAccessor.hpp"
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
  typedef
    CellAccessor<TGridGeometry, TMemory, D>
    SpecializedCellAccessor;

  typedef typename SpecializedCellAccessor::GridGeometry GridGeometry;
  typedef typename SpecializedCellAccessor::Memory       Memory;

  typedef Grid<SpecializedCellAccessor, D>              SpecializedGrid;
  typedef typename SpecializedGrid::VectorDi            VectorDi;
  typedef typename SpecializedGrid::CellAccessorFactory CellAccessorFactory;

  typedef SimulationParameters<Scalar, D> SpecializedSimulationParameters;

  typedef ParallelTopology<D> SpecializedParallelTopology;

  typedef
    Solvers::PoissonSolver<SpecializedCellAccessor, Scalar, D>
    PoissonSolver;

  typedef
    GhostCellsHandler<SpecializedGrid, D>
    SpecializedGhostCellsHandler;

public:
  MyTemplateSimulation() {}

  MyTemplateSimulation(MyTemplateSimulation const& other) = delete;

  ~MyTemplateSimulation() {}

  MyTemplateSimulation const&
  operator=(MyTemplateSimulation const& other) = delete;

  void
  initialize() {
    VectorDi localSize(_parallelTopology.localSize +
                       2 * VectorDi::Ones());

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
  }

  void
  iterate() {
    for (auto accessor : _grid.innerGrid) {
      typedef FghProcessing<SpecializedCellAccessor,
                            SpecializedSimulationParameters,
                            Scalar,
                            D> Fgh;
      Fgh::compute(accessor, _parameters, _dt);
    }

    _poissonSolver.solve();

    for (auto accessor : _grid.innerGrid) {
      typedef VelocityProcessing<SpecializedCellAccessor, Scalar, D> Velocity;
      Velocity::compute(accessor, _dt);
    }
  }

  Memory                          _memory;
  GridGeometry                    _gridGeometry;
  SpecializedGrid                 _grid;
  SpecializedSimulationParameters _parameters;
  SpecializedParallelTopology     _parallelTopology;
  Scalar                          _dt;
  PoissonSolver                   _poissonSolver;

  SpecializedGhostCellsHandler _ghostCellsHandler;
};
}
#endif
