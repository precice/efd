#ifndef FsiSimulation_CellAccessor_hpp
#define FsiSimulation_CellAccessor_hpp

#include <Uni/StructuredGrid/Basic/MultiIndex>

#include "StructuredMemory/Accessor.hpp"

namespace FsiSimulation {
template <typename TGridGeometry, typename TMemory, int D>
class CommonCellAccessor
  : public StructuredMemory::CommonAccessor<
      Uni::StructuredGrid::Basic::MultiIndex<D>, TMemory, D> {
public:
  typedef  StructuredMemory::CommonAccessor<
      Uni::StructuredGrid::Basic::MultiIndex<D>, TMemory, D> Base;
  typedef typename Base::VectorDi         VectorDi;
  typedef TMemory                         Memory;
  typedef typename Memory::Cell           Cell;
  typedef TGridGeometry                   GridGeometry;
  typedef typename GridGeometry::VectorDs VectorDs;

public:
  CommonCellAccessor(VectorDi const& i,
                     Memory*         memory,
                     GridGeometry*   gridGeometry)
    : Base(i, memory),
      _gridGeometry(gridGeometry) {}

  CommonCellAccessor(CommonCellAccessor const& other)
    : Base(other),
      _gridGeometry(other._gridGeometry) {}

  ~CommonCellAccessor() {}

  CommonCellAccessor&
  operator=(CommonCellAccessor const& other) {
    this->Base::operator=(other);
    _gridGeometry = other._gridGeometry;

    return *this;
  }

  VectorDs const&
  currentWidth() const {
    return _gridGeometry->cellWidth(this->indexValues());
  }

  VectorDs const&
  currentPosition() const {
    return _gridGeometry->cellPosition(this->indexValues());
  }

  VectorDs const&
  relativeWidth(VectorDi const& i) const {
    return _gridGeometry->cellWidth(this->indexValues() + i);
  }

  VectorDs const&
  relativePosition(VectorDi const& i) const {
    return _gridGeometry->cellPosition(this->indexValues() + i);
  }

  VectorDs const&
  absoluteWidth(VectorDi const& i) const {
    return _gridGeometry->cellWidth(i);
  }

  VectorDs const&
  absolutePosition(VectorDi const& i) const {
    return _gridGeometry->cellPosition(i);
  }

  VectorDs const&
  leftWidthInDimension(int const& dimension) const {
    return _gridGeometry->cellWidth(this->leftIndexInDimension(dimension));
  }

  VectorDs const&
  leftPositionInDimension(int const& dimension) const {
    return _gridGeometry->cellPosition(this->leftIndexInDimensiuon(dimension));
  }

  VectorDs
  rightWidthInDimension(int const& dimension) const {
    return _gridGeometry->cellWidth(this->rightIndexInDimension(dimension));
  }

  VectorDs
  rightPositionInDimension(int const& dimension) const {
    return _gridGeometry->cellPosition(this->rightIndexInDimensiuon(dimension));
  }

protected:
  GridGeometry* _gridGeometry;
};

template <typename TGridGeometry,
          typename TMemory,
          int D>
class CommonCellAccessor1D
  : public CommonCellAccessor<TGridGeometry, TMemory, D>,
    public StructuredMemory::CommonAccessor1D<
      typename CommonCellAccessor<TGridGeometry, TMemory, D>::Base, D> {
public:
  typedef  CommonCellAccessor<TGridGeometry, TMemory, D> Base;
  typedef typename Base::VectorDi                        VectorDi;
  typedef typename Base::Memory                          Memory;
  typedef typename Base::Cell                            Cell;
  typedef typename Base::GridGeometry                    GridGeometry;
  typedef typename Base::VectorDs                        VectorDs;

public:
  CommonCellAccessor1D(VectorDi const& i,
                       TMemory*        memory,
                       GridGeometry*   gridGeometry)
    : Base(i, memory, gridGeometry) {}

  CommonCellAccessor1D(CommonCellAccessor1D const& other)
    : Base(other) {}

  ~CommonCellAccessor1D() { this->Base::~Base(); }

  CommonCellAccessor1D&
  operator=(CommonCellAccessor1D const& other) {
    this->Base::operator=(other);

    return *this;
  }

  VectorDs
  leftWidth() const {
    return this->_gridGeometry->cellWidth(this->leftIndex());
  }

  VectorDs
  leftPosition() const {
    return this->_gridGeometry->cellPosition(this->leftIndex());
  }

  VectorDs
  rightWidth() const {
    return this->_gridGeometry->cellWidth(this->rightIndex());
  }

  VectorDs
  rightPosition() const {
    return this->_gridGeometry->cellPosition(this->rightIndex());
  }
};

template <typename TGridGeometry,
          typename TMemory,
          int D>
class CommonCellAccessor2D
  : public CommonCellAccessor1D<TGridGeometry, TMemory, D>,
    public StructuredMemory::CommonAccessor2D<
      typename CommonCellAccessor<TGridGeometry, TMemory, D>::Base, D> {
public:
  typedef  CommonCellAccessor1D<TGridGeometry, TMemory, D> Base;
  typedef typename Base::VectorDi                          VectorDi;
  typedef typename Base::Memory                            Memory;
  typedef typename Base::Cell                              Cell;
  typedef typename Base::GridGeometry                      GridGeometry;
  typedef typename Base::VectorDs                          VectorDs;

public:
  CommonCellAccessor2D(VectorDi const& i,
                       TMemory*        memory,
                       GridGeometry*   gridGeometry)
    : Base(i, memory, gridGeometry) {}

  CommonCellAccessor2D(CommonCellAccessor2D const& other)
    : Base(other) {}

  ~CommonCellAccessor2D() { this->Base::~Base(); }

  CommonCellAccessor2D&
  operator=(CommonCellAccessor2D const& other) {
    this->Base::operator=(other);

    return *this;
  }

  VectorDs
  bottomWidth() const {
    return this->_gridGeometry->cellWidth(this->bottomIndex());
  }

  VectorDs
  bottomPosition() const {
    return this->_gridGeometry->cellPosition(this->bottomIndex());
  }

  VectorDs
  topWidth() const {
    return this->_gridGeometry->cellWidth(this->topIndex());
  }

  VectorDs
  topPosition() const {
    return this->_gridGeometry->cellPosition(this->topIndex());
  }

  VectorDs
  leftBottomWidth() const {
    return this->_gridGeometry->cellWidth(this->leftBottomIndex());
  }

  VectorDs
  leftBottomPosition() const {
    return this->_gridGeometry->cellPosition(this->leftBottomIndex());
  }

  VectorDs
  leftTopWidth() const {
    return this->_gridGeometry->cellWidth(this->leftTopIndex());
  }

  VectorDs
  leftTopPosition() const {
    return this->_gridGeometry->cellPosition(this->leftTopIndex());
  }

  VectorDs
  rightBottomWidth() const {
    return this->_gridGeometry->cellWidth(this->rightBottomIndex());
  }

  VectorDs
  rightBottomPosition() const {
    return this->_gridGeometry->cellPosition(this->rightBottomIndex());
  }

  VectorDs
  rightTopWidth() const {
    return this->_gridGeometry->cellWidth(this->rightTopIndex());
  }

  VectorDs
  rightTopPosition() const {
    return this->_gridGeometry->cellPosition(this->rightTopIndex());
  }
};

template <typename TGridGeometry,
          typename TMemory,
          int D>
class CommonCellAccessor3D
  : public CommonCellAccessor2D<TGridGeometry, TMemory, D>,
    public StructuredMemory::CommonAccessor3D<
      typename CommonCellAccessor<TGridGeometry, TMemory, D>::Base, D> {
public:
  typedef  CommonCellAccessor2D<TGridGeometry, TMemory, D> Base;
  typedef typename Base::VectorDi                          VectorDi;
  typedef typename Base::Memory                            Memory;
  typedef typename Base::Cell                              Cell;
  typedef typename Base::GridGeometry                      GridGeometry;
  typedef typename Base::VectorDs                          VectorDs;

public:
  CommonCellAccessor3D(VectorDi const& i,
                       TMemory*        memory,
                       GridGeometry*   gridGeometry)
    : Base(i, memory, gridGeometry) {}

  CommonCellAccessor3D(CommonCellAccessor3D const& other) : Base(other) {}

  ~CommonCellAccessor3D() { this->Base::~Base(); }

  CommonCellAccessor3D&
  operator=(CommonCellAccessor3D const& other) {
    this->Base::operator=(other);

    return *this;
  }

  VectorDs
  backWidth() const {
    return this->_gridGeometry->cellWidth(this->backIndex());
  }

  VectorDs
  backPosition() const {
    return this->_gridGeometry->cellPosition(this->backIndex());
  }

  VectorDs
  frontWidth() const {
    return this->_gridGeometry->cellWidth(this->frontIndex());
  }

  VectorDs
  frontPosition() const {
    return this->_gridGeometry->cellPosition(this->frontIndex());
  }

  VectorDs
  leftBackWidth() const {
    return this->_gridGeometry->cellWidth(this->leftBackIndex());
  }

  VectorDs
  leftBackPosition() const {
    return this->_gridGeometry->cellPosition(this->leftBackIndex());
  }

  VectorDs
  leftFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->leftFrontIndex());
  }

  VectorDs
  leftFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->leftFrontIndex());
  }

  VectorDs
  rightBackWidth() const {
    return this->_gridGeometry->cellWidth(this->rightBackIndex());
  }

  VectorDs
  rightBackPosition() const {
    return this->_gridGeometry->cellPosition(this->rightBackIndex());
  }

  VectorDs
  rightFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->rightFrontIndex());
  }

  VectorDs
  rightFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->rightFrontIndex());
  }

  VectorDs
  bottomBackWidth() const {
    return this->_gridGeometry->cellWidth(this->bottomBackIndex());
  }

  VectorDs
  bottomBackPosition() const {
    return this->_gridGeometry->cellPosition(this->bottomBackIndex());
  }

  VectorDs
  bottomFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->bottomFrontIndex());
  }

  VectorDs
  bottomFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->bottomFrontIndex());
  }

  VectorDs
  topBackWidth() const {
    return this->_gridGeometry->cellWidth(this->topBackIndex());
  }

  VectorDs
  topBackPosition() const {
    return this->_gridGeometry->cellPosition(this->topBackIndex());
  }

  VectorDs
  topFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->topFrontIndex());
  }

  VectorDs
  topFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->topFrontIndex());
  }
};

template <typename TGridGeometry,
          typename TMemory,
          int D>
class CellAccessor
  : public CommonCellAccessor<TGridGeometry, TMemory, D> {
public:
  typedef  CommonCellAccessor<TGridGeometry, TMemory, D> Base;
  typedef typename Base::VectorDi                        VectorDi;
  typedef typename Base::Memory                          Memory;
  typedef typename Base::Cell                            Cell;
  typedef typename Base::GridGeometry                    GridGeometry;
  typedef typename Base::VectorDs                        VectorDs;

public:
  CellAccessor(VectorDi const& i,
               TMemory*        memory,
               GridGeometry*   gridGeometry)
    : Base(i, memory, gridGeometry) {}

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
          typename TMemory>
class CellAccessor<TGridGeometry, TMemory, 1>
  : public CommonCellAccessor1D<TGridGeometry, TMemory, 1> {
public:
  typedef  CommonCellAccessor1D<TGridGeometry, TMemory, 1> Base;
  typedef typename Base::VectorDi                          VectorDi;
  typedef typename Base::Memory                            Memory;
  typedef typename Base::Cell                              Cell;
  typedef typename Base::GridGeometry                      GridGeometry;
  typedef typename Base::VectorDs                          VectorDs;

public:
  CellAccessor(VectorDi const& i,
               TMemory*        memory,
               GridGeometry*   gridGeometry) : Base(i,
                                                    memory,
                                                    gridGeometry) {}

  CellAccessor(CellAccessor const& other) : Base(other) {}

  ~CellAccessor() { this->Base::~Base(); }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->Base::operator=(other);

    return *this;
  }
};

template <typename TGridGeometry, typename TMemory>
class CellAccessor<TGridGeometry, TMemory, 2>
  : public CommonCellAccessor2D<TGridGeometry, TMemory, 2> {
public:
  typedef  CommonCellAccessor2D<TGridGeometry, TMemory, 2> Base;
  typedef typename Base::VectorDi                          VectorDi;
  typedef typename Base::Memory                            Memory;
  typedef typename Base::Cell                              Cell;
  typedef typename Base::GridGeometry                      GridGeometry;
  typedef typename Base::VectorDs                          VectorDs;

public:
  CellAccessor(VectorDi const& i,
               TMemory*        memory,
               GridGeometry*   gridGeometry) : Base(i,
                                                    memory,
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
          typename TMemory>
class CellAccessor<TGridGeometry, TMemory, 3>
  : public CommonCellAccessor3D<TGridGeometry, TMemory, 3> {
public:
  typedef  CommonCellAccessor3D<TGridGeometry, TMemory, 3> Base;
  typedef typename Base::VectorDi                          VectorDi;
  typedef typename Base::Memory                            Memory;
  typedef typename Base::Cell                              Cell;
  typedef typename Base::GridGeometry                      GridGeometry;
  typedef typename Base::VectorDs                          VectorDs;

public:
  CellAccessor(VectorDi const& i,
               TMemory*        memory,
               GridGeometry*   gridGeometry) : Base(i,
                                                    memory,
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
