#ifndef FsiSimulation_MySimulation_hpp
#define FsiSimulation_MySimulation_hpp

#include "CellAccessor.hpp"
#include "Grid.hpp"
#include "Parameters.h"
#include "SimulationParameters.hpp"
#include "stencils/mystencils.hpp"

#include <Uni/Logging/macros>

namespace FsiSimulation {
class MySimulation {
public:
  virtual
  ~MySimulation() {}

  virtual void
  initialize(Parameters const& parameters) = 0;
};
template <typename GridGeometryType,
          typename CellHeapType,
          typename Scalar,
          int D>
class MyTemplateSimulation : public MySimulation {
public:
  typedef
    CellAccessor<GridGeometryType, CellHeapType, D>
    SpecializedCellAccessor;

  typedef typename SpecializedCellAccessor::GridGeometry GridGeometry;

  typedef typename SpecializedCellAccessor::CellHeap CellHeap;

  typedef Grid<SpecializedCellAccessor, D>              SpecializedGrid;
  typedef typename SpecializedGrid::VectorDi            VectorDi;
  typedef typename SpecializedGrid::CellAccessorFactory CellAccessorFactory;

  typedef SimulationParameters<Scalar, D> SpecializedSimulationParameters;

public:
  MyTemplateSimulation() {}

  MyTemplateSimulation(MyTemplateSimulation const& other) = delete;

  ~MyTemplateSimulation() {}

  MyTemplateSimulation const&
  operator=(MyTemplateSimulation const& other) = delete;

  void
  initialize(Parameters const& parameters) {
    VectorDi size(parameters.parallel.localSize);

    _cellHeap.allocate(size);

    typename GridGeometry::VectorDs cellSize;
    typename GridGeometry::VectorDi corner;
    cellSize(0) = parameters.geometry.lengthX / parameters.geometry.sizeX;
    cellSize(1) = parameters.geometry.lengthY / parameters.geometry.sizeY;
    corner(0)   = parameters.parallel.firstCorner[0];
    corner(1)   = parameters.parallel.firstCorner[1];

    if (D == 3) {
      cellSize(2) = parameters.geometry.lengthZ / parameters.geometry.sizeZ;
      corner(2)   = parameters.parallel.firstCorner[2];
    }

    _gridGeometry.cellSize(cellSize);
    _gridGeometry.corner(corner);

    CellAccessorFactory cellAccessorFactory(
      [&] (VectorDi const& i) {
        return SpecializedCellAccessor(i, &_cellHeap, &_gridGeometry);
      });

    _grid.initialize(size,
                     cellAccessorFactory);
    _parameters.re()    = parameters.flow.Re;
    _parameters.gamma() = parameters.solver.gamma;
    _parameters.g(0)    = parameters.environment.gx;
    _parameters.g(1)    = parameters.environment.gy;

    if (D == 3) {
      _parameters.g(2) = parameters.environment.gz;
    }
  }

  void
  iterate() {
    for (auto accessor : _grid) {
      typedef FghProcessing<SpecializedCellAccessor,
                            SpecializedSimulationParameters,
                            Scalar,
                            D> Fgh;
      Fgh::compute(accessor, _parameters, _dt);
    }
  }

private:
  CellHeap                        _cellHeap;
  GridGeometry                    _gridGeometry;
  SpecializedGrid                 _grid;
  SpecializedSimulationParameters _parameters;
  Scalar                          _dt;
};
}
#endif
