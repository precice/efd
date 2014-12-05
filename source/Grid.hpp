#ifndef FsiSimulation_Grid_hpp
#define FsiSimulation_Grid_hpp

#include <Uni/StructuredGrid/Basic/Grid>
#include <Uni/StructuredGrid/GridIterator>

#include <Uni/Logging/macros>

#include <array>

namespace FsiSimulation {
template <typename TCellAccessor, int D>
class Grid;
template <typename TCellAccessor, int D>
void
logGridInitializationInfo(Grid<TCellAccessor, D> const& grid);

template <typename TCellAccessor, int D>
class Grid : public Uni::StructuredGrid::Basic::Grid<TCellAccessor, D> {
public:
  typedef
    Uni::StructuredGrid::Basic::Grid<TCellAccessor, D>
    Base;

  typedef typename Base::Iterator   Iterator;
  typedef  typename Base::VectorDi  VectorDi;
  typedef typename Base::MultiIndex CellAccessor;
  typedef typename Base::Factory
    CellAccessorFactory;

public:
  Grid() {}

  Grid(Grid const& other) = delete;

  ~Grid() {}

  Grid&
  operator=(Grid const& other) = delete;

  void
  initialize(VectorDi const&            size,
             CellAccessorFactory const& factory) {
    this->Base::initialize(size,
                           VectorDi::Zero(),
                           VectorDi::Zero(),
                           factory);

    VectorDi indent(VectorDi::Ones());

    for (int i = 0; i < D; ++i) {
      VectorDi leftIndent(VectorDi::Zero());
      leftIndent(i) = 0;
      VectorDi rightIndent(VectorDi::Zero());
      rightIndent(i) = size(i) - 1;
      boundaries[i][0].initialize(size,
                                  leftIndent,
                                  rightIndent,
                                  factory);
      boundaries[i][1].initialize(size,
                                  rightIndent,
                                  leftIndent,
                                  factory);
      leftIndent     = indent;
      leftIndent(i)  = 0;
      rightIndent    = indent;
      rightIndent(i) = size(i) - rightIndent(i);
      indentedBoundaries[i][0].initialize(size,
                                          leftIndent,
                                          rightIndent,
                                          factory);
      indentedBoundaries[i][1].initialize(size,
                                          rightIndent,
                                          leftIndent,
                                          factory);
    }
    innerGrid.initialize(size, indent, indent, factory);
  }

public:
  Base                               innerGrid;
  std::array<std::array<Base, 2>, D> boundaries;
  std::array<std::array<Base, 2>, D> indentedBoundaries;

  friend void
  logGridInitializationInfo<TCellAccessor, D
                            >(Grid const& grid);
};

template <typename TCellAccessor, int D>
void
logGridInitializationInfo(Grid<TCellAccessor, D> const& grid) {
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

  for (int d = 0; d < D; ++d) {
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

#endif
