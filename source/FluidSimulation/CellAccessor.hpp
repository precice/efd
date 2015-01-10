#ifndef FsiSimulation_FluidSimulation_CellAccessor_hpp
#define FsiSimulation_FluidSimulation_CellAccessor_hpp

#include <Uni/StructuredGrid/Basic/Grid>

#include "StructuredMemory/Accessor.hpp"

#include <Uni/StructuredGrid/Basic/MultiIndex>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TDerived,
          typename TMemory,
          typename TGridGeometry,
          int TD>
class CommonCellAccessor
  : public StructuredMemory::CommonAccessor<
      Uni::StructuredGrid::Basic::MultiIndex<TD>, TMemory, TD> {
public:
  typedef StructuredMemory::CommonAccessor<
      Uni::StructuredGrid::Basic::MultiIndex<TD>, TMemory, TD> BaseType;
  typedef typename BaseType::VectorDi VectorDiType;

  typedef
    Uni::StructuredGrid::Basic::Grid<TDerived, TD>
    GridType;
  typedef TGridGeometry                       GridGeometryType;
  typedef typename GridGeometryType::VectorDs VectorDsType;
  typedef TMemory                             MemoryType;
  typedef typename MemoryType::Cell           CellType;

public:
  CommonCellAccessor(
    VectorDiType const& i,
    MemoryType*         memory,
    GridType const*     grid,
    GridGeometryType*   gridGeometry)
    : BaseType(i, memory),
      _grid(grid),
      _gridGeometry(gridGeometry) {}

  CommonCellAccessor(CommonCellAccessor const& other)
    : BaseType(other),
      _grid(other._grid),
      _gridGeometry(other._gridGeometry) {}

  ~CommonCellAccessor() {}

  CommonCellAccessor&
  operator=(CommonCellAccessor const& other) {
    this->BaseType::operator=(other);
    _grid         = other._grid;
    _gridGeometry = other._gridGeometry;

    return *this;
  }

  GridGeometryType const*
  gridGeometry() const {
    return _gridGeometry;
  }

  VectorDsType
  currentWidth() const {
    return _gridGeometry->cellWidth(this->indexValues()
                                    - this->_grid->leftIndent());
  }

  VectorDsType
  currentPosition() const {
    return _gridGeometry->cellPosition(this->indexValues()
                                       - this->_grid->leftIndent());
  }

  VectorDsType
  currentVelocityPosition(int const& dimension) const {
    VectorDsType result
      = _gridGeometry->cellPosition(this->indexValues()
                                    - this->_grid->leftIndent())
        + 0.5 * _gridGeometry->cellWidth(this->indexValues());
    result(dimension) +=
      0.5 * _gridGeometry->cellWidth(
        this->indexValues()) (dimension);
  }

  VectorDsType
  relativeWidth(VectorDiType const& i) const {
    return _gridGeometry->cellWidth(this->indexValues() + i
                                    - this->_grid->leftIndent());
  }

  VectorDsType
  relativePosition(VectorDiType const& i) const {
    return _gridGeometry->cellPosition(this->indexValues() + i
                                       - this->_grid->leftIndent());
  }

  VectorDsType
  relativeWidth(int const& dimension, int const& direction) const {
    return _gridGeometry->cellWidth(this->relativeIndex(dimension, direction)
                                    - this->_grid->leftIndent());
  }

  VectorDsType
  relativePosition(int const& dimension, int const& direction) const {
    return _gridGeometry->cellPosition(this->relativeIndex(dimension,
                                                           direction)
                                       - this->_grid->leftIndent());
  }

  VectorDsType
  absoluteWidth(VectorDiType const& i) const {
    return _gridGeometry->cellWidth(i
                                    - this->_grid->leftIndent());
  }

  VectorDsType
  absolutePosition(VectorDiType const& i) const {
    return _gridGeometry->cellPosition(i
                                       - this->_grid->leftIndent());
  }

  VectorDsType
  leftWidthInDimension(int const& dimension) const {
    return _gridGeometry->cellWidth(this->leftIndexInDimension(dimension)
                                    - this->_grid->leftIndent());
  }

  VectorDsType
  leftPositionInDimension(int const& dimension) const {
    return _gridGeometry->cellPosition(this->leftIndexInDimensiuon(dimension)
                                       - this->_grid->leftIndent());
  }

  VectorDsType
  rightWidthInDimension(int const& dimension) const {
    return _gridGeometry->cellWidth(this->rightIndexInDimension(dimension)
                                    - this->_grid->leftIndent());
  }

  VectorDsType
  rightPositionInDimension(int const& dimension) const {
    return _gridGeometry->cellPosition(this->rightIndexInDimensiuon(dimension)
                                       - this->_grid->leftIndent());
  }

protected:
  GridType const*   _grid;
  GridGeometryType* _gridGeometry;
};

template <typename TDerived,
          typename TMemory,
          typename TGridGeometry,
          int TD>
class CommonCellAccessor1D
  : public CommonCellAccessor<TDerived, TMemory, TGridGeometry, TD>,
    public StructuredMemory::CommonAccessor1D<
      typename CommonCellAccessor<TDerived, TMemory, TGridGeometry,
                                  TD>::BaseType,
      TD> {
public:
  typedef CommonCellAccessor<TDerived, TMemory, TGridGeometry, TD> BaseType;
  typedef typename BaseType::MemoryType                            MemoryType;
  typedef typename BaseType::GridType                              GridType;
  typedef typename BaseType::CellType                              CellType;
  typedef typename BaseType::GridGeometryType
    GridGeometryType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::VectorDsType VectorDsType;

public:
  CommonCellAccessor1D(VectorDiType const& i,
                       MemoryType*         memory,
                       GridType const*     grid,
                       GridGeometryType*   gridGeometry)
    : BaseType(i, memory, grid, gridGeometry) {}

  CommonCellAccessor1D(CommonCellAccessor1D const& other)
    : BaseType(other) {}

  ~CommonCellAccessor1D() {
    this->
    BaseType::~BaseType();
  }

  CommonCellAccessor1D&
  operator=(CommonCellAccessor1D const& other) {
    this->BaseType::operator=(other);

    return *this;
  }

  VectorDsType
  leftWidth() const {
    return this->_gridGeometry->cellWidth(this->leftIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  leftPosition() const {
    return this->_gridGeometry->cellPosition(this->leftIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  rightWidth() const {
    return this->_gridGeometry->cellWidth(this->rightIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  rightPosition() const {
    return this->_gridGeometry->cellPosition(this->rightIndex()
                                             - this->_grid->leftIndent());
  }
};

template <typename TDerived,
          typename TMemory,
          typename TGridGeometry,
          int TD>
class CommonCellAccessor2D
  : public CommonCellAccessor1D<TDerived, TMemory, TGridGeometry, TD>,
    public StructuredMemory::CommonAccessor2D<
      typename CommonCellAccessor<TDerived, TMemory, TGridGeometry,
                                  TD>::BaseType,
      TD> {
public:
  typedef CommonCellAccessor1D<TDerived, TMemory, TGridGeometry, TD> BaseType;
  typedef typename BaseType::MemoryType                              MemoryType;
  typedef typename BaseType::GridType                                GridType;
  typedef typename BaseType::CellType                                CellType;
  typedef typename BaseType::GridGeometryType
    GridGeometryType;
  typedef typename BaseType::VectorDiType
    VectorDiType;
  typedef typename BaseType::VectorDsType
    VectorDsType;

public:
  CommonCellAccessor2D(VectorDiType const& i,
                       MemoryType*         memory,
                       GridType const*     grid,
                       GridGeometryType*   gridGeometry)
    : BaseType(i, memory, grid, gridGeometry) {}

  CommonCellAccessor2D(CommonCellAccessor2D const& other)
    : BaseType(other) {}

  ~CommonCellAccessor2D() {
    this->
    BaseType::~BaseType();
  }

  CommonCellAccessor2D&
  operator=(CommonCellAccessor2D const& other) {
    this->BaseType::operator=(other);

    return *this;
  }

  VectorDsType
  bottomWidth() const {
    return this->_gridGeometry->cellWidth(this->bottomIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  bottomPosition() const {
    return this->_gridGeometry->cellPosition(this->bottomIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  topWidth() const {
    return this->_gridGeometry->cellWidth(this->topIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  topPosition() const {
    return this->_gridGeometry->cellPosition(this->topIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  leftBottomWidth() const {
    return this->_gridGeometry->cellWidth(this->leftBottomIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  leftBottomPosition() const {
    return this->_gridGeometry->cellPosition(this->leftBottomIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  leftTopWidth() const {
    return this->_gridGeometry->cellWidth(this->leftTopIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  leftTopPosition() const {
    return this->_gridGeometry->cellPosition(this->leftTopIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  rightBottomWidth() const {
    return this->_gridGeometry->cellWidth(this->rightBottomIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  rightBottomPosition() const {
    return this->_gridGeometry->cellPosition(this->rightBottomIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  rightTopWidth() const {
    return this->_gridGeometry->cellWidth(this->rightTopIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  rightTopPosition() const {
    return this->_gridGeometry->cellPosition(this->rightTopIndex()
                                             - this->_grid->leftIndent());
  }
};

template <typename TDerived,
          typename TMemory,
          typename TGridGeometry,
          int TD>
class CommonCellAccessor3D
  : public CommonCellAccessor2D<TDerived, TMemory, TGridGeometry, TD>,
    public StructuredMemory::CommonAccessor3D<
      typename CommonCellAccessor<TDerived, TMemory, TGridGeometry,
                                  TD>::BaseType,
      TD> {
public:
  typedef CommonCellAccessor2D<TDerived, TMemory, TGridGeometry, TD> BaseType;
  typedef typename BaseType::MemoryType                              MemoryType;
  typedef typename BaseType::GridType                                GridType;
  typedef typename BaseType::GridGeometryType
    GridGeometryType;
  typedef typename BaseType::CellType CellType;
  typedef typename BaseType::VectorDiType
    VectorDiType;
  typedef typename BaseType::VectorDsType
    VectorDsType;

public:
  CommonCellAccessor3D(VectorDiType const& i,
                       MemoryType*         memory,
                       GridType const*     grid,
                       GridGeometryType*   gridGeometry)
    : BaseType(i, memory, grid, gridGeometry) {}

  CommonCellAccessor3D(CommonCellAccessor3D const& other) : BaseType(other) {}

  ~CommonCellAccessor3D() {
    this->
    BaseType::~BaseType();
  }

  CommonCellAccessor3D&
  operator=(CommonCellAccessor3D const& other) {
    this->BaseType::operator=(other);

    return *this;
  }

  VectorDsType
  backWidth() const {
    return this->_gridGeometry->cellWidth(this->backIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  backPosition() const {
    return this->_gridGeometry->cellPosition(this->backIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  frontWidth() const {
    return this->_gridGeometry->cellWidth(this->frontIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  frontPosition() const {
    return this->_gridGeometry->cellPosition(this->frontIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  leftBackWidth() const {
    return this->_gridGeometry->cellWidth(this->leftBackIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  leftBackPosition() const {
    return this->_gridGeometry->cellPosition(this->leftBackIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  leftFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->leftFrontIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  leftFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->leftFrontIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  rightBackWidth() const {
    return this->_gridGeometry->cellWidth(this->rightBackIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  rightBackPosition() const {
    return this->_gridGeometry->cellPosition(this->rightBackIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  rightFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->rightFrontIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  rightFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->rightFrontIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  bottomBackWidth() const {
    return this->_gridGeometry->cellWidth(this->bottomBackIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  bottomBackPosition() const {
    return this->_gridGeometry->cellPosition(this->bottomBackIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  bottomFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->bottomFrontIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  bottomFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->bottomFrontIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  topBackWidth() const {
    return this->_gridGeometry->cellWidth(this->topBackIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  topBackPosition() const {
    return this->_gridGeometry->cellPosition(this->topBackIndex()
                                             - this->_grid->leftIndent());
  }

  VectorDsType
  topFrontWidth() const {
    return this->_gridGeometry->cellWidth(this->topFrontIndex()
                                          - this->_grid->leftIndent());
  }

  VectorDsType
  topFrontPosition() const {
    return this->_gridGeometry->cellPosition(this->topFrontIndex()
                                             - this->_grid->leftIndent());
  }
};

template <typename TMemory,
          typename TGridGeometry,
          int TD>
class CellAccessor
  : public CommonCellAccessor<
      CellAccessor<TMemory, TGridGeometry, TD>,
      TMemory, TGridGeometry, TD> {
public:
  typedef CommonCellAccessor<CellAccessor, TMemory, TGridGeometry, TD> BaseType;
  typedef typename BaseType::VectorDiType
    VectorDiType;
  typedef typename BaseType::MemoryType
    MemoryType;
  typedef typename BaseType::GridType GridType;
  typedef typename BaseType::CellType CellType;
  typedef typename BaseType::GridGeometryType
    GridGeometryType;

public:
  CellAccessor(VectorDiType const& i,
               MemoryType*         memory,
               GridType const*     grid,
               GridGeometryType*   gridGeometry)
    : BaseType(i, memory, grid, gridGeometry) {}

  CellAccessor(CellAccessor const& other)
    : BaseType(other) {}

  ~CellAccessor() {
    this->
    BaseType::~BaseType();
  }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->BaseType::operator=(other);

    return *this;
  }
};

template <typename TMemory,
          typename TGridGeometry>
class CellAccessor<TMemory, TGridGeometry, 1>
  : public CommonCellAccessor1D<
      CellAccessor<TMemory, TGridGeometry, 1>,
      TMemory, TGridGeometry, 1> {
public:
  typedef CommonCellAccessor1D<
      CellAccessor<TMemory, TGridGeometry, 1>,
      TMemory, TGridGeometry, 1>
    BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType   MemoryType;
  typedef typename BaseType::GridType     GridType;
  typedef typename BaseType::CellType     CellType;
  typedef typename BaseType::GridGeometryType
    GridGeometryType;

public:
  CellAccessor(VectorDiType const& i,
               MemoryType*         memory,
               GridType const*     grid,
               GridGeometryType*   gridGeometry)
    : BaseType(i, memory, grid, gridGeometry) {}

  CellAccessor(CellAccessor const& other) : BaseType(other) {}

  ~CellAccessor() {
    this->
    BaseType::~BaseType();
  }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->BaseType::operator=(other);

    return *this;
  }
};

template <typename TMemory,
          typename TGridGeometry>
class CellAccessor<TMemory, TGridGeometry, 2>
  : public CommonCellAccessor2D<
      CellAccessor<TMemory, TGridGeometry, 2>,
      TMemory, TGridGeometry, 2> {
public:
  typedef CommonCellAccessor2D<
      CellAccessor<TMemory, TGridGeometry, 2>,
      TMemory, TGridGeometry, 2>
    BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType   MemoryType;
  typedef typename BaseType::GridType     GridType;
  typedef typename BaseType::CellType     CellType;
  typedef typename BaseType::GridGeometryType
    GridGeometryType;

public:
  CellAccessor(VectorDiType const& i,
               MemoryType*         memory,
               GridType const*     grid,
               GridGeometryType*   gridGeometry)
    : BaseType(i, memory, grid, gridGeometry) {}

  CellAccessor(CellAccessor const& other) : BaseType(other) {}

  ~CellAccessor() {
    this->
    BaseType::~BaseType();
  }

  CellAccessor&
  operator=(CellAccessor const& other) {
    this->BaseType::operator=(other);

    return *this;
  }
};

template <typename TMemory,
          typename TGridGeometry>
class CellAccessor<TMemory, TGridGeometry, 3>
  : public CommonCellAccessor3D<
      CellAccessor<TMemory, TGridGeometry, 3>,
      TMemory, TGridGeometry, 3> {
public:
  typedef CommonCellAccessor3D<
      CellAccessor<TMemory, TGridGeometry, 3>,
      TMemory, TGridGeometry, 3>
    BaseType;
  typedef typename BaseType::VectorDiType VectorDiType;
  typedef typename BaseType::MemoryType   MemoryType;
  typedef typename BaseType::GridType     GridType;
  typedef typename BaseType::CellType     CellType;
  typedef typename BaseType::GridGeometryType
    GridGeometryType;

public:
  CellAccessor(VectorDiType const& i,
               MemoryType*         memory,
               GridType const*     grid,
               GridGeometryType*   gridGeometry)
    : BaseType(i, memory, grid, gridGeometry) {}

  CellAccessor(CellAccessor const& other) : BaseType(other) {}

  ~CellAccessor() {
    this->
    BaseType::~BaseType();
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
