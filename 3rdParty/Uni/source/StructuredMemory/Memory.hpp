#pragma once

#include "Accessor.hpp"
#include "Pointers.hpp"

#include <Uni/StructuredGrid/Basic/Grid>
#include <Uni/StructuredGrid/Basic/MultiIndex>

#include <Eigen/Core>

namespace Uni {
namespace StructuredMemory {
template <typename TCell, int TD>
class Memory {
public:
  typedef TCell                     Cell;
  typedef Eigen::Matrix<int, TD, 1> VectorDi;
  typedef Pointers<Cell, TD>        SpecializedPointers;

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

  VectorDi const&
  size() const {
    return _pointers.size();
  }

  Cell*
  data() const {
    return _pointers.data();
  }

  VectorDi const&
  indexShift() const {
    return _indexShift;
  }

  void
  indexShift(VectorDi const index) {
    _indexShift = index;
  }

  Cell*
  getCell(VectorDi const& index) {
    return _pointers.getCell(index);
  }

  Cell*
  getCell(VectorDi const& index,
          VectorDi const& indexShift) {
    return _pointers.getCell(index + indexShift);
  }

  Cell*
  getCellWithShiftIndex(VectorDi const& index) {
    return _pointers.getCell(index + _indexShift);
  }

  Cell*
  operator()(VectorDi const& index) {
    return getCell(index);
  }

  Cell*
  operator()(VectorDi const& index,
             VectorDi const& indexShift) {
    return getCell(index + indexShift);
  }

private:
  SpecializedPointers _pointers;
  VectorDi            _indexShift;
};

template <typename TCell, int TD>
class IterableMemory : public Memory<TCell, TD> {
public:
  typedef Memory<TCell, TD>                  Base;
  typedef typename Base::Cell                Cell;
  typedef typename Base::VectorDi            VectorDi;
  typedef typename Base::SpecializedPointers SpecializedPointers;

  struct MultiIndexTraits {
    using Type = Uni::StructuredGrid::Basic::MultiIndex<MultiIndexTraits>;

    enum {
      Dimensions = TD
    };
  };

  typedef
    Accessor<Uni::StructuredGrid::Basic::MultiIndex<MultiIndexTraits>, Base, TD>
    SpecializedAccessor;

  typedef
    Uni::StructuredGrid::Basic::Grid<SpecializedAccessor>
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
