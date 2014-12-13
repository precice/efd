#ifndef FsiSimulation_CellAccessor_hpp
#define FsiSimulation_CellAccessor_hpp

#include <Uni/StructuredGrid/Basic/MultiIndex>

#include "StructuredMemory/Accessor.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template<typename TGridGeometry, typename TMemory, int D>
class CommonCellAccessor
: public StructuredMemory::CommonAccessor<
Uni::StructuredGrid::Basic::MultiIndex<D>, TMemory, D> {
public:
  typedef StructuredMemory::CommonAccessor<
  Uni::StructuredGrid::Basic::MultiIndex<D>, TMemory, D> BaseType;
  typedef typename BaseType::VectorDi VectorDiType;
  typedef TMemory MemoryType;
  typedef typename MemoryType::Cell CellType;
  typedef TGridGeometry GridGeometryType;
  typedef typename GridGeometryType::VectorDs VectorDsType;

public:
  CommonCellAccessor(VectorDiType const& i,
  MemoryType* memory,
  GridGeometryType* gridGeometry)
  : BaseType(i, memory),
    _gridGeometry(gridGeometry) {
  }

  CommonCellAccessor(CommonCellAccessor const& other)
  : BaseType(other),
    _gridGeometry(other._gridGeometry) {
  }

  ~CommonCellAccessor() {
  }

  CommonCellAccessor&
  operator=(CommonCellAccessor const& other) {
    this->BaseType::operator=(other);
    _gridGeometry = other._gridGeometry;

    return *this;
  }

  VectorDsType
  currentWidth() const {
    return _gridGeometry->cellWidth(this->indexValues());
  }

  VectorDsType
  currentPosition() const {
    return _gridGeometry->cellPosition(this->indexValues());
  }

  VectorDsType
  relativeWidth(VectorDiType const& i) const {
    return _gridGeometry->cellWidth(this->indexValues() + i);
  }

  VectorDsType
  relativePosition(VectorDiType const& i) const {
    return _gridGeometry->cellPosition(this->indexValues() + i);
  }

  VectorDsType
  absoluteWidth(VectorDiType const& i) const {
    return _gridGeometry->cellWidth(i);
  }

  VectorDsType
  absolutePosition(VectorDiType const& i) const {
    return _gridGeometry->cellPosition(i);
  }

  VectorDsType
  leftWidthInDimension(int const& dimension) const {
    return _gridGeometry->cellWidth(this->leftIndexInDimension(dimension));
  }

  VectorDsType
  leftPositionInDimension(int const& dimension) const {
    return _gridGeometry->cellPosition(this->leftIndexInDimensiuon(dimension));
  }

  VectorDsType
  rightWidthInDimension(int const& dimension) const {
    return _gridGeometry->cellWidth(this->rightIndexInDimension(dimension));
  }

  VectorDsType
  rightPositionInDimension(int const& dimension) const {
    return _gridGeometry->cellPosition(this->rightIndexInDimensiuon(dimension));
  }

protected:
  GridGeometryType* _gridGeometry;
};

template<typename TGridGeometry,
typename TMemory,
int D>
class CommonCellAccessor1D
: public CommonCellAccessor<TGridGeometry, TMemory, D>,
  public StructuredMemory::CommonAccessor1D<
  typename CommonCellAccessor<TGridGeometry, TMemory, D>::BaseType, D> {
public:
  typedef CommonCellAccessor<TGridGeometry, TMemory, D> BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType MemoryType;
  typedef typename BaseType::CellType CellType;
  typedef typename BaseType::GridGeometryType GridGeometryType;
  typedef typename BaseType::VectorDsType VectorDsType;

public:
  CommonCellAccessor1D(VectorDiType const& i,
  TMemory* memory,
  GridGeometryType* gridGeometry)
  : BaseType(i, memory, gridGeometry) {
  }

  CommonCellAccessor1D(CommonCellAccessor1D const& other)
  : BaseType(other) {
  }

  ~CommonCellAccessor1D() {
    this->BaseType::~BaseType();
  }

  CommonCellAccessor1D&
  operator=(CommonCellAccessor1D const& other) {
    this->BaseType::operator=(other);

    return *this;
  }

  VectorDsType
  leftWidth() const {
    return this->_gridGeometry->cellWidth(this->leftIndex());
  }

  VectorDsType
  leftPosition() const {
    return this->_gridGeometry->cellPosition(this->leftIndex());
  }

  VectorDsType
  rightWidth() const {
    return this->_gridGeometry->cellWidth(this->rightIndex());
  }

  VectorDsType
  rightPosition() const {
    return this->_gridGeometry->cellPosition(this->rightIndex());
  }
};

template<typename TGridGeometry,
typename TMemory,
int D>
class CommonCellAccessor2D
: public CommonCellAccessor1D<TGridGeometry, TMemory, D>,
  public StructuredMemory::CommonAccessor2D<
  typename CommonCellAccessor<TGridGeometry, TMemory, D>::BaseType, D> {
public:
  typedef CommonCellAccessor1D<TGridGeometry, TMemory, D> BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType MemoryType;
  typedef typename BaseType::CellType CellType;
  typedef typename BaseType::GridGeometryType GridGeometryType;
  typedef typename BaseType::VectorDsType VectorDsType;

public:
  CommonCellAccessor2D(VectorDiType const& i,
  TMemory* memory,
  GridGeometryType* gridGeometry)
  : BaseType(i, memory, gridGeometry) {
  }

  CommonCellAccessor2D(CommonCellAccessor2D const& other)
  : BaseType(other) {
  }

  ~CommonCellAccessor2D() {
    this->BaseType::~BaseType();
  }

  CommonCellAccessor2D&
  operator=(CommonCellAccessor2D const& other) {
    this->BaseType::operator=(other);

    return *this;
  }

  VectorDsType
  bottomWidth() const {
    return this->_gridGeometry->cellWidth(this->bottomIndex());
  }

  VectorDsType
  bottomPosition() const {
    return this->_gridGeometry->cellPosition(this->bottomIndex());
  }

  VectorDsType
  topWidth() const {
    return this->_gridGeometry->cellWidth(this->topIndex());
  }

  VectorDsType
  topPosition() const {
    return this->_gridGeometry->cellPosition(this->topIndex());
  }

  VectorDsType
  leftBottomWidth() const {
    return this->_gridGeometry->cellWidth(this->leftBottomIndex());
  }

  VectorDsType
  leftBottomPosition() const {
    return this->_gridGeometry->cellPosition(this->leftBottomIndex());
  }

  VectorDsType
  leftTopWidth() const {
    return this->_gridGeometry->cellWidth(this->leftTopIndex());
  }

  VectorDsType
  leftTopPosition() const {
    return this->_gridGeometry->cellPosition(this->leftTopIndex());
  }

  VectorDsType
  rightBottomWidth() const {
    return this->_gridGeometry->cellWidth(this->rightBottomIndex());
  }

  VectorDsType
  rightBottomPosition() const {
    return this->_gridGeometry->cellPosition(this->rightBottomIndex());
  }

  VectorDsType
  rightTopWidth() const {
    return this->_gridGeometry->cellWidth(this->rightTopIndex());
  }

  VectorDsType
  rightTopPosition() const {
    return this->_gridGeometry->cellPosition(this->rightTopIndex());
  }
};

template<typename TGridGeometry,
typename TMemory,
int D>
class CommonCellAccessor3D
: public CommonCellAccessor2D<TGridGeometry, TMemory, D>,
  public StructuredMemory::CommonAccessor3D<
  typename CommonCellAccessor<TGridGeometry, TMemory, D>::BaseType, D> {
public:
  typedef CommonCellAccessor2D<TGridGeometry, TMemory, D> BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType MemoryType;
  typedef typename BaseType::CellType CellType;
  typedef typename BaseType::GridGeometryType GridGeometryType;
  typedef typename BaseType::VectorDsType VectorDsType;

public:
  CommonCellAccessor3D(VectorDiType const& i,
  TMemory* memory,
  GridGeometryType* gridGeometry)
  : BaseType(i, memory, gridGeometry) {
  }

  CommonCellAccessor3D(CommonCellAccessor3D const& other) : BaseType(other) {
  }

  ~CommonCellAccessor3D() {
    this->BaseType::~BaseType();
  }

  CommonCellAccessor3D&
  operator=(CommonCellAccessor3D const& other) {
    this->BaseType::operator=(other);

    return *this;
  }

  VectorDsType
  backWidth() const {
    return this->_gridGeometry->cellWidth(this->backIndex());
  }

  VectorDsType
  backPosition() const {
    return this->_gridGeometry->cellPosition(this->backIndex());
  }

  VectorDsType
  frontWidth() const {
    return this->_gridGeometry->cellWidth(this->frontIndex());
  }

  VectorDsType
  frontPosition() const {
    return this->_gridGeometry->cellPosition(this->frontIndex());
  }

  VectorDsType
  leftBackWidth() const {
    return this->_gridGeometry->cellWidth(this->leftBackIndex());
  }

  VectorDsType
  leftBackPosition() const {
    return this->_gridGeometry->cellPosition(this->leftBackIndex());
  }

  VectorDsType
  leftFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->leftFrontIndex());
  }

  VectorDsType
  leftFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->leftFrontIndex());
  }

  VectorDsType
  rightBackWidth() const {
    return this->_gridGeometry->cellWidth(this->rightBackIndex());
  }

  VectorDsType
  rightBackPosition() const {
    return this->_gridGeometry->cellPosition(this->rightBackIndex());
  }

  VectorDsType
  rightFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->rightFrontIndex());
  }

  VectorDsType
  rightFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->rightFrontIndex());
  }

  VectorDsType
  bottomBackWidth() const {
    return this->_gridGeometry->cellWidth(this->bottomBackIndex());
  }

  VectorDsType
  bottomBackPosition() const {
    return this->_gridGeometry->cellPosition(this->bottomBackIndex());
  }

  VectorDsType
  bottomFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->bottomFrontIndex());
  }

  VectorDsType
  bottomFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->bottomFrontIndex());
  }

  VectorDsType
  topBackWidth() const {
    return this->_gridGeometry->cellWidth(this->topBackIndex());
  }

  VectorDsType
  topBackPosition() const {
    return this->_gridGeometry->cellPosition(this->topBackIndex());
  }

  VectorDsType
  topFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->topFrontIndex());
  }

  VectorDsType
  topFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->topFrontIndex());
  }
};

template<typename TGridGeometry,
typename TMemory,
int D>
class CellAccessor
: public CommonCellAccessor<TGridGeometry, TMemory, D> {
public:
  typedef CommonCellAccessor<TGridGeometry, TMemory, D> BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType MemoryType;
  typedef typename BaseType::CellType CellType;
  typedef typename BaseType::GridGeometryType GridGeometryType;

public:
  CellAccessor(VectorDiType const& i,
  TMemory* memory,
  GridGeometryType* gridGeometry)
  : BaseType(i, memory, gridGeometry) {
  }

  CellAccessor(CellAccessor const& other)
  : BaseType(other) {
  }

  ~CellAccessor() {
    this->BaseType::~BaseType();
  }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->BaseType::operator=(other);

    return *this;
  }
};

template<typename TGridGeometry,
typename TMemory>
class CellAccessor<TGridGeometry, TMemory, 1>
: public CommonCellAccessor1D<TGridGeometry, TMemory, 1> {
public:
  typedef CommonCellAccessor1D<TGridGeometry, TMemory, 1> BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType MemoryType;
  typedef typename BaseType::CellType CellType;
  typedef typename BaseType::GridGeometryType GridGeometryType;

public:
  CellAccessor(VectorDiType const& i,
  TMemory* memory,
  GridGeometryType* gridGeometry) : BaseType(i,
  memory,
  gridGeometry) {
  }

  CellAccessor(CellAccessor const& other) : BaseType(other) {
  }

  ~CellAccessor() {
    this->BaseType::~BaseType();
  }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->BaseType::operator=(other);

    return *this;
  }
};

template<typename TGridGeometry, typename TMemory>
class CellAccessor<TGridGeometry, TMemory, 2>
: public CommonCellAccessor2D<TGridGeometry, TMemory, 2> {
public:
  typedef CommonCellAccessor2D<TGridGeometry, TMemory, 2> BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType MemoryType;
  typedef typename BaseType::CellType CellType;
  typedef typename BaseType::GridGeometryType GridGeometryType;

public:
  CellAccessor(VectorDiType const& i,
  TMemory* memory,
  GridGeometryType* gridGeometry) : BaseType(i,
  memory,
  gridGeometry) {
  }

  CellAccessor(CellAccessor const& other) : BaseType(other) {
  }

  ~CellAccessor() {
    this->BaseType::~BaseType();
  }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->BaseType::operator=(other);

    return *this;
  }
};

template<typename TGridGeometry,
typename TMemory>
class CellAccessor<TGridGeometry, TMemory, 3>
: public CommonCellAccessor3D<TGridGeometry, TMemory, 3> {
public:
  typedef CommonCellAccessor3D<TGridGeometry, TMemory, 3> BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType MemoryType;
  typedef typename BaseType::CellType CellType;
  typedef typename BaseType::GridGeometryType GridGeometryType;

public:
  CellAccessor(VectorDiType const& i,
  TMemory* memory,
  GridGeometryType* gridGeometry)
  : BaseType(i, memory, gridGeometry) {
  }

  CellAccessor(CellAccessor const& other) : BaseType(other) {
  }

  ~CellAccessor() {
    this->BaseType::~BaseType();
  }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->BaseType::operator=(other);

    return *this;
  }
};
}
}
#endif
