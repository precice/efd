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

  using FactoryType = std::function
                      <CellAccessorType(BaseType const*, VectorDiType const&)>;

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

    typename BaseType::Factory cell_accessor_factory(
      [&] (VectorDiType const& index) {
        return _factory(this, index);
      });

    this->BaseType::initialize(size,
                               VectorDiType::Zero(),
                               VectorDiType::Zero(),
                               cell_accessor_factory);

    VectorDiType indent(VectorDiType::Ones());

    for (int i = 0; i < Dimensions; ++i) {
      VectorDiType leftIndent(VectorDiType::Zero());
      leftIndent(i) = 0;
      VectorDiType rightIndent(VectorDiType::Zero());
      rightIndent(i) = size(i) - 1;

      typename BaseType::Factory cell_accessor_factory_1(
        [&] (VectorDiType const& index) {
          return _factory(&boundaries[i][0], index);
        });
      boundaries[i][0].initialize(size,
                                  leftIndent,
                                  rightIndent,
                                  cell_accessor_factory_1);

      typename BaseType::Factory cell_accessor_factory_2(
        [&] (VectorDiType const& index) {
          return _factory(&boundaries[i][1], index);
        });
      boundaries[i][1].initialize(size,
                                  rightIndent,
                                  leftIndent,
                                  cell_accessor_factory_2);
      leftIndent     = indent;
      leftIndent(i)  = 0;
      rightIndent    = indent;
      rightIndent(i) = size(i) - rightIndent(i);

      typename BaseType::Factory cell_accessor_factory_3(
        [&] (VectorDiType const& index) {
          return _factory(&indentedBoundaries[i][1], index);
        });
      indentedBoundaries[i][0].initialize(size,
                                          leftIndent,
                                          rightIndent,
                                          cell_accessor_factory_3);

      typename BaseType::Factory cell_accessor_factory_4(
        [&] (VectorDiType const& index) {
          return _factory(&indentedBoundaries[i][1], index);
        });
      indentedBoundaries[i][1].initialize(size,
                                          rightIndent,
                                          leftIndent,
                                          cell_accessor_factory_4);
    }

    typename BaseType::Factory cell_accessor_factory_5(
      [&] (VectorDiType const& index) {
        return _factory(&innerGrid, index);
      });
    innerGrid.initialize(size, indent, indent, cell_accessor_factory_5);
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
