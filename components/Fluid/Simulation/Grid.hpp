#pragma once

#include <Uni/Logging/macros>
#include <Uni/StructuredGrid/Basic/Grid>

#include <array>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class Grid;

template <typename TSolverTraits>
void
logGridInitializationInfo(Grid<TSolverTraits> const& grid);

template <typename TSolverTraits>
class Grid
  : public Uni::StructuredGrid::Basic::Grid
    <typename TSolverTraits::CellAccessorType> {
public:
  using BaseType = Uni::StructuredGrid::Basic::Grid
                   <typename TSolverTraits::CellAccessorType>;

  using CellAccessorType =  typename TSolverTraits::CellAccessorType;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

  using VectorDiType = typename BaseType::VectorDi;

  using FactoryType = typename BaseType::Factory;

  using IteratorType = typename BaseType::Iterator;

public:
  Grid() {}

  Grid(Grid const& other) = delete;

  ~Grid() {}

  Grid&
  operator=(Grid const& other) = delete;

  void
  initialize(VectorDiType const& size,
             FactoryType const&  factory) {
    _factory = factory;

    this->BaseType::initialize(size,
                               VectorDiType::Zero(),
                               VectorDiType::Zero(),
                               _factory);

    VectorDiType indent(VectorDiType::Ones());

    for (int i = 0; i < Dimensions; ++i) {
      VectorDiType leftIndent(VectorDiType::Zero());
      leftIndent(i) = 0;
      VectorDiType rightIndent(VectorDiType::Zero());
      rightIndent(i) = size(i) - 1;

      boundaries[i][0].initialize(size,
                                  leftIndent,
                                  rightIndent,
                                  _factory);

      boundaries[i][1].initialize(size,
                                  rightIndent,
                                  leftIndent,
                                  _factory);
      leftIndent     = indent;
      leftIndent(i)  = 0;
      rightIndent    = indent;
      rightIndent(i) = size(i) - rightIndent(i);

      indentedBoundaries[i][0].initialize(size,
                                          leftIndent,
                                          rightIndent,
                                          _factory);

      indentedBoundaries[i][1].initialize(size,
                                          rightIndent,
                                          leftIndent,
                                          _factory);
    }

    innerGrid.initialize(size, indent, indent, _factory);
  }

public:
  BaseType                                        innerGrid;
  FactoryType                                     _factory;
  std::array<std::array<BaseType, 2>, Dimensions> boundaries;
  std::array<std::array<BaseType, 2>, Dimensions> indentedBoundaries;

  friend void
  logGridInitializationInfo<TSolverTraits>(Grid const& grid);
};

template <typename TSolverTraits>
void
logGridInitializationInfo(Grid<TSolverTraits> const& grid) {
  enum {
    Dimensions = TSolverTraits::Dimensions
  };

  INFO << "Grid size: "
       << grid.size().transpose() << "\n"
       << "Grid left indent: "
       << grid.leftIndent().transpose() << "\n"
       << "Grid right indent: "
       << grid.rightIndent().transpose() << "\n"
       << "Inner grid size: "
       << grid.innerGrid.size().transpose() << "\n"
       << "Inner grid left indent: "
       << grid.innerGrid.leftIndent().transpose() << "\n"
       << "Inner grid right indent: "
       << grid.innerGrid.rightIndent().transpose() << "\n";

  for (int d = 0; d < Dimensions; ++d) {
    INFO
      << 2 * d << " grid size: "
      << grid.boundaries[d][0].size().transpose() << "\n"
      << 2 * d << " grid left indent: "
      << grid.boundaries[d][0].leftIndent().transpose() << "\n"
      << 2 * d << " grid right indent: "
      << grid.boundaries[d][0].rightIndent().transpose() << "\n"
      << 2 * d + 1 << " grid size: "
      << grid.boundaries[d][1].size().transpose() << "\n"
      << 2 * d + 1 << " grid left indent: "
      << grid.boundaries[d][1].leftIndent().transpose() << "\n"
      << 2 * d + 1 << " grid right indent: "
      << grid.boundaries[d][1].rightIndent().transpose() << "\n";
  }
}
}
}
