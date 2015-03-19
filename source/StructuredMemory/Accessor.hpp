#ifndef FsiSimulation_StructuredMemory_Accessor_hpp
#define FsiSimulation_StructuredMemory_Accessor_hpp

namespace FsiSimulation {
namespace StructuredMemory {
template <typename TDerived, int TD>
class CommonAccessor1D;

template <typename TDerived, int TD>
class CommonAccessor2D;

template <typename TDerived, int TD>
class CommonAccessor3D;

template <typename TMultiIndex, typename TMemory, int TD>
class CommonAccessor : public TMultiIndex {
  friend class CommonAccessor1D<CommonAccessor<TMultiIndex, TMemory, TD>, TD>;

  friend class CommonAccessor2D<CommonAccessor<TMultiIndex, TMemory, TD>, TD>;

  friend class CommonAccessor3D<CommonAccessor<TMultiIndex, TMemory, TD>, TD>;

public:
  typedef TMultiIndex             Base;
  typedef typename Base::VectorDi VectorDi;
  typedef TMemory                 Memory;
  typedef typename Memory::Cell   Cell;

public:
  CommonAccessor(VectorDi const& i,
                 Memory*         memory) : Base(i),
    _memory(memory) {}

  CommonAccessor(CommonAccessor const& other) : Base(other),
    _memory(other._memory) {}

  ~CommonAccessor() {}

  CommonAccessor&
  operator=(CommonAccessor const& other) {
    this->Base::operator=(other);
    _memory = other._memory;

    return *this;
  }

  VectorDi
  currentIndex() const {
    VectorDi result(this->indexValues());

    return result;
  }

  Cell*
  currentCell() const {
    return (*_memory)(this->indexValues());
  }

  VectorDi
  relativeIndex(VectorDi const& i) const {
    VectorDi result(this->indexValues() + i);

    return result;
  }

  Cell*
  relativeCell(VectorDi const& i) const {
    return (*_memory)(this->indexValues() + i);
  }

  VectorDi
  relativeIndex(int const& dimension, int const& direction,
                int const& distance = 1) const {
    VectorDi result(this->indexValues());

    if (direction == 0) {
      result(dimension) -= distance;
    } else {
      result(dimension) += distance;
    }

    return result;
  }

  Cell*
  relativeCell(int const& dimension, int const& direction,
               int const& distance = 1) const {
    return (*_memory)(relativeIndex(dimension, direction, distance));
  }

  Cell*
  absoluteCell(VectorDi const& i) const {
    return (*_memory)(i);
  }

  VectorDi
  leftIndexInDimension(int const& dimension) const {
    VectorDi result(this->indexValues());
    result(dimension) -= 1;

    return result;
  }

  Cell*
  leftCellInDimension(int const& dimension) const {
    return (*_memory)(leftIndexInDimension(dimension));
  }

  VectorDi
  rightIndexInDimension(int const& dimension) const {
    VectorDi result(this->indexValues());
    result(dimension) += 1;

    return result;
  }

  Cell*
  rightCellInDimension(int const& dimension) const {
    return (*_memory)(rightIndexInDimension(dimension));
  }

  VectorDi
  leftRightIndexInDimensions(int const& dimension,
                             int const& dimension2) const {
    VectorDi result(this->indexValues());
    result(dimension)  -= 1;
    result(dimension2) += 1;

    return result;
  }

  Cell*
  leftRightCellInDimensions(int const& dimension,
                            int const& dimension2) const {
    return (*_memory)(leftRightIndexInDimensions(dimension,
                                                 dimension2));
  }

protected:
  Memory* _memory;
};

template <typename TDerived, typename T>
inline TDerived&
cast(T* pointer) {
  return *reinterpret_cast<TDerived*>(pointer);
}

template <typename TDerived, typename T>
inline TDerived const&
cast(T const* pointer) {
  return *reinterpret_cast<TDerived const*>(pointer);
}

template <typename TDerived, int TD>
class CommonAccessor1D {
public:
  typedef TDerived                   Derived;
  typedef typename Derived::VectorDi VectorDi;
  typedef typename Derived::Memory   Memory;
  typedef typename Derived::Cell     Cell;

public:
  VectorDi
  leftIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) -= 1;

    return result;
  }

  Cell*
  leftCell() const {
    return (*cast<Derived>(this)._memory)(leftIndex());
  }

  VectorDi
  rightIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) += 1;

    return result;
  }

  Cell*
  rightCell() const {
    return (*cast<Derived>(this)._memory)(rightIndex());
  }
};

template <typename TDerived, int TD>
class CommonAccessor2D {
public:
  typedef TDerived                   Derived;
  typedef typename Derived::VectorDi VectorDi;
  typedef typename Derived::Memory   Memory;
  typedef typename Derived::Cell     Cell;

public:
  VectorDi
  bottomIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(1) -= 1;

    return result;
  }

  Cell*
  bottomCell() const {
    return (*cast<Derived>(this)._memory)(bottomIndex());
  }

  VectorDi
  topIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(1) += 1;

    return result;
  }

  Cell*
  topCell() const {
    return (*cast<Derived>(this)._memory)(topIndex());
  }

  VectorDi
  leftBottomIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) -= 1;
    result(1) -= 1;

    return result;
  }

  Cell*
  leftBottomCell() const {
    return (*cast<Derived>(this)._memory)(leftBottomIndex());
  }

  VectorDi
  leftTopIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) -= 1;
    result(1) += 1;

    return result;
  }

  Cell*
  leftTopCell() const {
    return (*cast<Derived>(this)._memory)(leftTopIndex());
  }

  VectorDi
  rightBottomIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) += 1;
    result(1) -= 1;

    return result;
  }

  Cell*
  rightBottomCell() const {
    return (*cast<Derived>(this)._memory)(rightBottomIndex());
  }

  VectorDi
  rightTopIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) += 1;
    result(1) += 1;

    return result;
  }

  Cell*
  rightTopCell() const {
    return (*cast<Derived>(this)._memory)(rightTopIndex());
  }
};

template <typename TDerived, int TD>
class CommonAccessor3D {
public:
  typedef TDerived                   Derived;
  typedef typename Derived::VectorDi VectorDi;
  typedef typename Derived::Memory   Memory;
  typedef typename Derived::Cell     Cell;

public:
  VectorDi
  backIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(2) -= 1;

    return result;
  }

  Cell*
  backCell() const {
    return (*cast<Derived>(this)._memory)(backIndex());
  }

  VectorDi
  frontIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(2) -= 1;

    return result;
  }

  Cell*
  frontCell() const {
    return (*cast<Derived>(this)._memory)(frontIndex());
  }

  VectorDi
  leftBackIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) -= 1;
    result(2) -= 1;

    return result;
  }

  Cell*
  leftBackCell() const {
    return (*cast<Derived>(this)._memory)(leftBackIndex());
  }

  VectorDi
  leftFrontIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) -= 1;
    result(2) += 1;

    return result;
  }

  Cell*
  leftFrontCell() const {
    return (*cast<Derived>(this)._memory)(leftFrontIndex());
  }

  VectorDi
  rightBackIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) += 1;
    result(2) -= 1;

    return result;
  }

  Cell*
  rightBackCell() const {
    return (*cast<Derived>(this)._memory)(
      rightBackIndex());
  }

  VectorDi
  rightFrontIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(0) += 1;
    result(2) += 1;

    return result;
  }

  Cell*
  rightFrontCell() const {
    return (*cast<Derived>(this)._memory)(rightFrontIndex());
  }

  VectorDi
  bottomBackIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(1) -= 1;
    result(2) -= 1;

    return result;
  }

  Cell*
  bottomBackCell() const {
    return (*cast<Derived>(this)._memory)(bottomBackIndex());
  }

  VectorDi
  bottomFrontIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(1) -= 1;
    result(2) += 1;

    return result;
  }

  Cell*
  bottomFrontCell() const {
    return (*cast<Derived>(this)._memory)(
      bottomFrontIndex());
  }

  VectorDi
  topBackIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(1) += 1;
    result(2) -= 1;

    return result;
  }

  Cell*
  topBackCell() const {
    return (*cast<Derived>(this)._memory)(topBackIndex());
  }

  VectorDi
  topFrontIndex() const {
    VectorDi result(cast<Derived>(this).indexValues());
    result(1) += 1;
    result(2) += 1;

    return result;
  }

  Cell*
  topFrontCell() const {
    return (*cast<Derived>(this)._memory)(topFrontIndex());
  }
};

template <typename TMultiIndex, typename TMemory, int TD>
class Accessor
  : public CommonAccessor<TMultiIndex, TMemory, TD> {
public:
  typedef CommonAccessor<TMultiIndex, TMemory, TD> Base;
  typedef typename Base::VectorDi                  VectorDi;
  typedef typename Base::Memory                    Memory;
  typedef typename Base::Cell                      Cell;

public:
  Accessor(VectorDi const& i,
           TMemory*        memory) : Base(i, memory) {}

  Accessor(Accessor const& other) : Base(other) {}

  ~Accessor() {
    this->
    Base::~Base();
  }

  Accessor&
  operator=(Accessor const& other) {
    this->Base::operator=(other);

    return *this;
  }
};

template <typename TMultiIndex, typename TMemory>
class Accessor<TMultiIndex, TMemory, 1>
  : public CommonAccessor<TMultiIndex, TMemory, 1>,
    public CommonAccessor1D<CommonAccessor<TMultiIndex, TMemory, 1>, 1> {
public:

  enum {
    Dimensions = 1
  };

  typedef CommonAccessor<TMultiIndex, TMemory, 1> Base;
  typedef typename Base::VectorDi                 VectorDi;
  typedef typename Base::Memory                   Memory;
  typedef typename Base::Cell                     Cell;

public:
  Accessor(VectorDi const& i,
           TMemory*        memory) : Base(i, memory) {}

  Accessor(Accessor const& other) : Base(other) {}

  ~Accessor() {
    this->
    Base::~Base();
  }

  Accessor&
  operator=(Accessor const& other) {
    this->Base::operator=(other);

    return *this;
  }
};

template <typename TMultiIndex, typename TMemory>
class Accessor<TMultiIndex, TMemory, 2>
  : public CommonAccessor<TMultiIndex, TMemory, 2>,
    public CommonAccessor1D<CommonAccessor<TMultiIndex, TMemory, 2>, 2>,
    public CommonAccessor2D<CommonAccessor<TMultiIndex, TMemory, 2>, 2> {
public:

  enum {
    Dimensions = 2
  };

  typedef CommonAccessor<TMultiIndex, TMemory, 2> Base;
  typedef typename Base::VectorDi                 VectorDi;
  typedef typename Base::Memory                   Memory;
  typedef typename Base::Cell                     Cell;

public:
  Accessor(VectorDi const& i,
           TMemory*        memory) : Base(i, memory) {}

  Accessor(Accessor const& other) : Base(other) {}

  ~Accessor() {
    this->
    Base::~Base();
  }

  Accessor&
  operator=(Accessor const& other) {
    this->Base::operator=(other);

    return *this;
  }
};

template <typename TMultiIndex, typename TMemory>
class Accessor<TMultiIndex, TMemory, 3>
  : public CommonAccessor<TMultiIndex, TMemory, 3>,
    public CommonAccessor1D<CommonAccessor<TMultiIndex, TMemory, 3>, 3>,
    public CommonAccessor2D<CommonAccessor<TMultiIndex, TMemory, 3>, 3>,
    public CommonAccessor3D<CommonAccessor<TMultiIndex, TMemory, 3>, 3> {
public:

  enum {
    Dimensions = 3
  };

  typedef CommonAccessor<TMultiIndex, TMemory, 3> Base;

  typedef typename Base::VectorDi VectorDi;
  typedef typename Base::Memory   Memory;
  typedef typename Base::Cell     Cell;

public:
  Accessor(VectorDi const& i,
           TMemory*        memory) : Base(i, memory) {}

  Accessor(Accessor const& other) : Base(other) {}

  ~Accessor() {
    this->
    Base::~Base();
  }

  Accessor&
  operator=(Accessor const& other) {
    this->Base::operator=(other);

    return *this;
  }
};
}
}
#endif
