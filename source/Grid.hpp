#ifndef FsiSimulation_Grid_hpp
#define FsiSimulation_Grid_hpp

#include <Uni/StructuredGrid/Basic/Grid>
#include <Uni/StructuredGrid/GridIterator>

namespace FsiSimulation {
template <typename CellAccessor, int D>
class Grid : public Uni::StructuredGrid::Basic::Grid<CellAccessor, D> {
public:
  typedef Uni::StructuredGrid::Basic::Grid<CellAccessor, D> Base;
  typedef typename Base::Iterator                           Iterator;
  typedef  typename Base::VectorDi                          VectorDi;
  typedef typename Base::Factory                            CellAccessorFactory;

public:
  Grid() {}

  Grid(Grid const& other) = delete;

  ~Grid() {}

  Grid&
  operator=(Grid const& other) = delete;

  void
  initialize(VectorDi const& size,
             CellAccessorFactory const&  factory) {
    VectorDi indent;
    indent.fill(1);
    this->Base::initialize(size,
                           indent,
                           indent,
                           factory);

    for (int i = 0; i < D; ++i) {
      auto leftIndent(indent);
      leftIndent(i) = 0;
      auto rightIndent(indent);
      rightIndent(i) = size(i) - rightIndent(i);
      boundaries[2 * i].initialize(size,
                                   leftIndent,
                                   rightIndent,
                                   factory);
      boundaries[2 * i + 1].initialize(size,
                                       rightIndent,
                                       leftIndent,
                                       factory);
    }
  }

public:
  Base boundaries[2 * D];
};
}

#endif
