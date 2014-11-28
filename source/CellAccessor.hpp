#ifndef FsiSimulation_CellAccessor_hpp
#define FsiSimulation_CellAccessor_hpp

#include "CellHeap.hpp"

#include <Uni/StructuredGrid/Basic/MultiIndex>

namespace FsiSimulation {
template <typename TGridGeometry,
          typename TCellHeap,
          int D>
class CommonCellAccessor : public Uni::StructuredGrid::Basic::MultiIndex<D> {
public:
  typedef Uni::StructuredGrid::Basic::MultiIndex<D> Base;
  typedef typename Base::VectorDi                   VectorDi;
  typedef TCellHeap                                 CellHeap;
  typedef typename CellHeap::Cell                   Cell;
  typedef TGridGeometry                             GridGeometry;
  typedef typename GridGeometry::VectorDs           VectorDs;

public:
  CommonCellAccessor(VectorDi const& i,
                     TCellHeap*      cellHeap,
                     GridGeometry*   gridGeometry)
    : Base(i),
      _cellHeap(cellHeap),
      _gridGeometry(gridGeometry) {}

  CommonCellAccessor(CommonCellAccessor const& other)
    : Base(other),
      _cellHeap(other._cellHeap),
      _gridGeometry(
        other._gridGeometry) {}

  ~CommonCellAccessor() {}

  CommonCellAccessor&
  operator=(CommonCellAccessor const& other) {
    this->Base::operator=(other);
    _cellHeap     = other._cellHeap;
    _gridGeometry = other._gridGeometry;

    return *this;
  }

  Cell*
  currentCell() const {
    return (* _cellHeap)(Base::indexValues());
  }

  VectorDs
  currentSize() const {
    return _gridGeometry->cellSize(Base::indexValues());
  }

  VectorDs
  currentPosition() const {
    return _gridGeometry->cellPosition(Base::indexValues());
  }

  Cell*
  relativeCell(VectorDi const& i) const {
    return (* _cellHeap)(Base::indexValues() + i);
  }

  Cell*
  relativeSize(VectorDi const& i) const {
    return _gridGeometry(Base::indexValues() + i);
  }

  Cell*
  relativePosition(VectorDi const& i) const {
    return _gridGeometry(Base::indexValues() + i);
  }

protected:
  TCellHeap*    _cellHeap;
  GridGeometry* _gridGeometry;
};

template <typename TGridGeometry,
          typename TCellHeap,
          int D>
class CommonCellAccessor1D
  : public CommonCellAccessor<TGridGeometry, TCellHeap, D> {
public:
  typedef  CommonCellAccessor<TGridGeometry, TCellHeap, D> Base;
  typedef typename Base::VectorDi                          VectorDi;
  typedef typename Base::CellHeap                          CellHeap;
  typedef typename Base::Cell                              Cell;
  typedef typename Base::GridGeometry                      GridGeometry;
  typedef typename Base::VectorDs                          VectorDs;

public:
  CommonCellAccessor1D(VectorDi const& i,
                       TCellHeap*      cellHeap,
                       GridGeometry*   gridGeometry)
    : Base(i, cellHeap, gridGeometry) {}

  CommonCellAccessor1D(CommonCellAccessor1D const& other)
    : Base(other) {}

  ~CommonCellAccessor1D() { this->Base::~Base(); }

  CommonCellAccessor1D&
  operator=(CommonCellAccessor1D const& other) {
    this->Base::operator=(other);

    return *this;
  }

  VectorDi
  leftIndex() const {
    VectorDi result(Base::indexValues());
    result(0) -= 1;

    return result;
  }

  Cell*
  leftCell() const {
    return (*this->_cellHeap)(leftIndex());
  }

  VectorDs
  leftSize() const { return this->_gridGeometry->cellSize(leftIndex()); }

  VectorDs
  leftPosition() const {
    return this->_gridGeometry->cellPosition(leftIndex());
  }

  VectorDi
  rightIndex() const {
    VectorDi result(Base::indexValues());
    result(0) += 1;

    return result;
  }

  Cell*
  rightCell() const { return (*this->_cellHeap)(rightIndex()); }

  VectorDs
  rightSize() const {
    return this->_gridGeometry->cellSize(rightIndex());
  }

  VectorDs
  rightPosition() const {
    return this->_gridGeometry->cellPosition(rightIndex());
  }
};

template <typename TGridGeometry,
          typename TCellHeap,
          int D>
class CommonCellAccessor2D
  : public CommonCellAccessor1D<TGridGeometry, TCellHeap, D> {
public:
  typedef  CommonCellAccessor1D<TGridGeometry, TCellHeap, D> Base;
  typedef typename Base::VectorDi                            VectorDi;
  typedef typename Base::CellHeap                            CellHeap;
  typedef typename Base::Cell                                Cell;
  typedef typename Base::GridGeometry                        GridGeometry;
  typedef typename Base::VectorDs                            VectorDs;

public:
  CommonCellAccessor2D(VectorDi const& i,
                       TCellHeap*      cellHeap,
                       GridGeometry*   gridGeometry)
    : Base(i, cellHeap, gridGeometry) {}

  CommonCellAccessor2D(CommonCellAccessor2D const& other)
    : Base(other) {}

  ~CommonCellAccessor2D() { this->Base::~Base(); }

  CommonCellAccessor2D&
  operator=(CommonCellAccessor2D const& other) {
    this->Base::operator=(other);

    return *this;
  }

  VectorDi
  bottomIndex() const {
    VectorDi result(Base::indexValues());
    result(1) -= 1;

    return result;
  }

  Cell*
  bottomCell() const { return (*this->_cellHeap)(bottomIndex()); }

  VectorDs
  bottomSize() const {
    return this->_gridGeometry->cellSize(bottomIndex());
  }

  VectorDs
  bottomPosition() const {
    return this->_gridGeometry->cellPosition(bottomIndex());
  }

  VectorDi
  topIndex() const {
    VectorDi result(Base::indexValues());
    result(1) += 1;

    return result;
  }

  Cell*
  topCell() const {
    return (*this->_cellHeap)(topIndex());
  }

  VectorDs
  topSize() const {
    return this->_gridGeometry->cellSize(topIndex());
  }

  VectorDs
  topPosition() const {
    return this->_gridGeometry->cellPosition(topIndex());
  }

  VectorDi
  leftBottomIndex() const {
    VectorDi result(Base::indexValues());
    result(0) -= 1;
    result(1) -= 1;

    return result;
  }

  Cell*
  leftBottomCell() const { return (*this->_cellHeap)(leftBottomIndex()); }

  VectorDs
  leftBottomSize() const {
    return this->_gridGeometry->cellSize(leftBottomIndex());
  }

  VectorDs
  leftBottomPosition() const {
    return this->_gridGeometry->cellPosition(leftBottomIndex());
  }

  VectorDi
  leftTopIndex() const {
    VectorDi result(Base::indexValues());
    result(0) -= 1;
    result(1) += 1;

    return result;
  }

  Cell*
  leftTopCell() const { return (*this->_cellHeap)(leftTopIndex()); }

  VectorDs
  leftTopSize() const {
    return this->_gridGeometry->cellSize(leftTopIndex());
  }

  VectorDs
  leftTopPosition() const {
    return this->_gridGeometry->cellPosition(leftTopIndex());
  }

  VectorDi
  rightBottomIndex() const {
    VectorDi result(Base::indexValues());
    result(0) += 1;
    result(1) -= 1;

    return result;
  }

  Cell*
  rightBottomCell() const { return (*this->_cellHeap)(rightBottomIndex()); }

  VectorDs
  rightBottomSize() const {
    return this->_gridGeometry->cellSize(rightBottomIndex());
  }

  VectorDs
  rightBottomPosition() const {
    return this->_gridGeometry->cellPosition(rightBottomIndex());
  }

  VectorDi
  rightTopIndex() const {
    VectorDi result(Base::indexValues());
    result(0) += 1;
    result(1) += 1;

    return result;
  }

  Cell*
  rightTopCell() const { return (*this->_cellHeap)(rightTopIndex()); }

  VectorDs
  rightTopSize() const {
    return this->_gridGeometry->cellSize(rightTopIndex());
  }

  VectorDs
  rightTopPosition() const {
    return this->_gridGeometry->cellPosition(rightTopIndex());
  }
};

template <typename TGridGeometry,
          typename TCellHeap,
          int D>
class CommonCellAccessor3D
  : public CommonCellAccessor2D<TGridGeometry, TCellHeap, D> {
public:
  typedef  CommonCellAccessor2D<TGridGeometry, TCellHeap, D> Base;
  typedef typename Base::VectorDi                            VectorDi;
  typedef typename Base::CellHeap                            CellHeap;
  typedef typename Base::Cell                                Cell;
  typedef typename Base::GridGeometry                        GridGeometry;
  typedef typename Base::VectorDs                            VectorDs;

public:
  CommonCellAccessor3D(VectorDi const& i,
                       TCellHeap*      cellHeap,
                       GridGeometry*   gridGeometry)
    : Base(i, cellHeap, gridGeometry) {}

  CommonCellAccessor3D(CommonCellAccessor3D const& other) : Base(other) {}

  ~CommonCellAccessor3D() { this->Base::~Base(); }

  CommonCellAccessor3D&
  operator=(CommonCellAccessor3D const& other) {
    this->Base::operator=(other);

    return *this;
  }

  VectorDi
  backIndex() const {
    VectorDi result(Base::indexValues());
    result(2) -= 1;

    return result;
  }

  Cell*
  backCell() const {
    return (*this->_cellHeap)(backIndex());
  }

  VectorDs
  backSize() const {
    return this->_gridGeometry->cellSize(backIndex());
  }

  VectorDs
  backPosition() const {
    return this->_gridGeometry->cellPosition(backIndex());
  }

  VectorDi
  frontIndex() const {
    VectorDi result(Base::indexValues());
    result(2) -= 1;

    return result;
  }

  Cell*
  frontCell() const {
    return (*this->_cellHeap)(frontIndex());
  }

  VectorDs
  frontSize() const {
    return this->_gridGeometry->cellSize(frontIndex());
  }

  VectorDs
  frontPosition() const {
    return this->_gridGeometry->cellPosition(frontIndex());
  }

  VectorDi
  leftBackIndex() const {
    VectorDi result(Base::indexValues());
    result(0) -= 1;
    result(2) -= 1;

    return result;
  }

  Cell*
  leftBackCell() const { return (*this->_cellHeap)(leftBackIndex()); }

  VectorDs
  leftBackSize() const {
    return this->_gridGeometry->cellSize(leftBackIndex());
  }

  VectorDs
  leftBackPosition() const {
    return this->_gridGeometry->cellPosition(leftBackIndex());
  }

  VectorDi
  leftFrontIndex() const {
    VectorDi result(Base::indexValues());
    result(0) -= 1;
    result(2) += 1;

    return result;
  }

  Cell*
  leftFrontCell() const { return (*this->_cellHeap)(leftFrontIndex()); }

  VectorDs
  leftFrontSize() const {
    return this->_gridGeometry->cellSize(leftFrontIndex());
  }

  VectorDs
  leftFrontPosition() const {
    return this->_gridGeometry->cellPosition(leftFrontIndex());
  }

  VectorDi
  rightBackIndex() const {
    VectorDi result(Base::indexValues());
    result(0) += 1;
    result(2) -= 1;

    return result;
  }

  Cell*
  rightBackCell() const { return (*this->_cellHeap)(rightBackIndex()); }

  VectorDs
  rightBackSize() const {
    return this->_gridGeometry->cellSize(rightBackIndex());
  }

  VectorDs
  rightBackPosition() const {
    return this->_gridGeometry->cellPosition(rightBackIndex());
  }

  VectorDi
  rightFrontIndex() const {
    VectorDi result(Base::indexValues());
    result(0) += 1;
    result(2) += 1;

    return result;
  }

  Cell*
  rightFrontCell() const { return (*this->_cellHeap)(rightFrontIndex()); }

  VectorDs
  rightFrontSize() const {
    return this->_gridGeometry->cellSize(rightFrontIndex());
  }

  VectorDs
  rightFrontPosition() const {
    return this->_gridGeometry->cellPosition(rightFrontIndex());
  }

  VectorDi
  bottomBackIndex() const {
    VectorDi result(Base::indexValues());
    result(1) -= 1;
    result(2) -= 1;

    return result;
  }

  Cell*
  bottomBackCell() const { return (*this->_cellHeap)(bottomBackIndex()); }

  VectorDs
  bottomBackSize() const {
    return this->_gridGeometry->cellSize(bottomBackIndex());
  }

  VectorDs
  bottomBackPosition() const {
    return this->_gridGeometry->cellPosition(bottomBackIndex());
  }

  VectorDi
  bottomFrontIndex() const {
    VectorDi result(Base::indexValues());
    result(1) -= 1;
    result(2) += 1;

    return result;
  }

  Cell*
  bottomFrontCell() const { return (*this->_cellHeap)(bottomFrontIndex()); }

  VectorDs
  bottomFrontSize() const {
    return this->_gridGeometry->cellSize(bottomFrontIndex());
  }

  VectorDs
  bottomFrontPosition() const {
    return this->_gridGeometry->cellPosition(bottomFrontIndex());
  }

  VectorDi
  topBackIndex() const {
    VectorDi result(Base::indexValues());
    result(1) += 1;
    result(2) -= 1;

    return result;
  }

  Cell*
  topBackCell() const { return (*this->_cellHeap)(topBackIndex()); }

  VectorDs
  topBackSize() const {
    return this->_gridGeometry->cellSize(topBackIndex());
  }

  VectorDs
  topBackPosition() const {
    return this->_gridGeometry->cellPosition(topBackIndex());
  }

  VectorDi
  topFrontIndex() const {
    VectorDi result(Base::indexValues());
    result(1) += 1;
    result(2) += 1;

    return result;
  }

  Cell*
  topFrontCell() const { return (*this->_cellHeap)(topFrontIndex()); }

  VectorDs
  topFrontSize() const {
    return this->_gridGeometry->cellSize(topFrontIndex());
  }

  VectorDs
  topFrontPosition() const {
    return this->_gridGeometry->cellPosition(topFrontIndex());
  }
};

template <typename TGridGeometry,
          typename TCellHeap,
          int D>
class CellAccessor
  : public CommonCellAccessor<TGridGeometry, TCellHeap, D> {
public:
  typedef  CommonCellAccessor<TGridGeometry, TCellHeap, D> Base;
  typedef typename Base::VectorDi                          VectorDi;
  typedef typename Base::CellHeap                          CellHeap;
  typedef typename Base::Cell                              Cell;
  typedef typename Base::GridGeometry                      GridGeometry;
  typedef typename Base::VectorDs                          VectorDs;

public:
  CellAccessor(VectorDi const& i,
               TCellHeap*      cellHeap,
               GridGeometry*   gridGeometry)
    : Base(i, cellHeap, gridGeometry) {}

  CellAccessor(CellAccessor const& other)
    : Base(other) {}

  ~CellAccessor() { this->Base::~Base(); }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->Base::operator=(other);

    return *this;
  }
};

template <typename TGridGeometry,
          typename TCellHeap>
class CellAccessor<TGridGeometry, TCellHeap, 1>
  : public CommonCellAccessor1D<TGridGeometry, TCellHeap, 1> {
public:
  typedef  CommonCellAccessor1D<TGridGeometry, TCellHeap, 1> Base;
  typedef typename Base::VectorDi                            VectorDi;
  typedef typename Base::CellHeap                            CellHeap;
  typedef typename Base::Cell                                Cell;
  typedef typename Base::GridGeometry                        GridGeometry;
  typedef typename Base::VectorDs                            VectorDs;

public:
  CellAccessor(VectorDi const& i,
               TCellHeap*      cellHeap,
               GridGeometry*   gridGeometry) : Base(i,
                                                    cellHeap,
                                                    gridGeometry) {}

  CellAccessor(CellAccessor const& other) : Base(other) {}

  ~CellAccessor() { this->Base::~Base(); }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->Base::operator=(other);

    return *this;
  }
};

template <typename TGridGeometry,
          typename TCellHeap>
class CellAccessor<TGridGeometry, TCellHeap, 2>
  : public CommonCellAccessor2D<TGridGeometry, TCellHeap, 2> {
public:
  typedef  CommonCellAccessor2D<TGridGeometry, TCellHeap, 2> Base;
  typedef typename Base::VectorDi                            VectorDi;
  typedef typename Base::CellHeap                            CellHeap;
  typedef typename Base::Cell                                Cell;
  typedef typename Base::GridGeometry                        GridGeometry;
  typedef typename Base::VectorDs                            VectorDs;

public:
  CellAccessor(VectorDi const& i,
               TCellHeap*      cellHeap,
               GridGeometry*   gridGeometry) : Base(i,
                                                    cellHeap,
                                                    gridGeometry) {}

  CellAccessor(CellAccessor const& other) : Base(other) {}

  ~CellAccessor() { this->Base::~Base(); }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->Base::operator=(other);

    return *this;
  }
};

template <typename TGridGeometry,
          typename TCellHeap>
class CellAccessor<TGridGeometry, TCellHeap, 3>
  : public CommonCellAccessor3D<TGridGeometry, TCellHeap, 3> {
public:
  typedef  CommonCellAccessor3D<TGridGeometry, TCellHeap, 3> Base;
  typedef typename Base::VectorDi                            VectorDi;
  typedef typename Base::CellHeap                            CellHeap;
  typedef typename Base::Cell                                Cell;
  typedef typename Base::GridGeometry                        GridGeometry;
  typedef typename Base::VectorDs                            VectorDs;

public:
  CellAccessor(VectorDi const& i,
               TCellHeap*      cellHeap,
               GridGeometry*   gridGeometry) : Base(i,
                                                    cellHeap,
                                                    gridGeometry) {}

  CellAccessor(CellAccessor const& other) : Base(other) {}

  ~CellAccessor() { this->Base::~Base(); }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->Base::operator=(other);

    return *this;
  }
};
}
#endif
