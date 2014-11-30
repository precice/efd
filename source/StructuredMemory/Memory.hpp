#ifndef FsiSimulation_StructuredMemory_Memory_hpp
#define FsiSimulation_StructuredMemory_Memory_hpp

#include "Accessor.hpp"
#include "Pointers.hpp"

#include <Uni/StructuredGrid/Basic/Grid>
#include <Uni/StructuredGrid/Basic/MultiIndex>

#include <Eigen/Core>

namespace FsiSimulation {
namespace StructuredMemory {
template <typename TCell, int D>
class Memory {
public:
  typedef TCell                    Cell;
  typedef Eigen::Matrix<int, D, 1> VectorDi;
  typedef Pointers<Cell, D>        SpecializedPointers;

public:
  Memory() {}

  Memory(Memory const& other) = delete;

  ~Memory() {}

  Memory&
  operator=(Memory const& other) = delete;

  void
  allocate(VectorDi const& size) {
    _pointers.allocate(size);
  }

  void
  release() {
    _pointers.release();
  }

  Cell*
  getCell(VectorDi const& index) {
    return _pointers.getCell(index);
  }

  Cell*
  operator()(VectorDi const& index) {
    return getCell(index);
  }

private:
  SpecializedPointers _pointers;
};

template <typename TCell, int D>
class IterableMemory : public Memory<TCell, D> {
public:
  typedef Memory<TCell, D>                   Base;
  typedef typename Base::Cell                Cell;
  typedef typename Base::VectorDi            VectorDi;
  typedef typename Base::SpecializedPointers SpecializedPointers;

  typedef
    Accessor<Uni::StructuredGrid::Basic::MultiIndex<D>, Base, D>
    SpecializedAccessor;

  typedef
    Uni::StructuredGrid::Basic::Grid<SpecializedAccessor, D>
    Grid;
  typedef typename Grid::Iterator Iterator;
  typedef typename Grid::Factory  AccessorFactory;

public:
  IterableMemory() : Base() {}

  IterableMemory(IterableMemory const& other) = delete;

  ~IterableMemory() {}

  IterableMemory&
  operator=(IterableMemory const& other) = delete;

  Iterator
  begin() {
    return _grid.begin();
  }

  Iterator
  end() {
    return _grid.end();
  }

  void
  allocate(VectorDi const& size) {
    AccessorFactory factory(
      [&] (VectorDi const& i) {
        return SpecializedAccessor(i, this);
      });
    _grid.initialize(size,
                     VectorDi::Zero(),
                     VectorDi::Zero(),
                     factory);
    this->Base::allocate(size);
  }

private:
  Grid _grid;
};
}
}
#endif
