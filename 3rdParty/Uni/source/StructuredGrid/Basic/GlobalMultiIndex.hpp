#ifndef Uni_StructuredGrid_Basic_GlobalMultiIndex_hpp
#define Uni_StructuredGrid_Basic_GlobalMultiIndex_hpp

#include <Uni/StructuredGrid/MultiIndex>

#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>

#include <array>

#define Uni_Structured_Grid_Basic_COMPUTE_GLOBAL_INDEX_2D(x, y, z, xSize, ySize) \
  (x) + (y) * (xSize) + (z) * (xSize) * (ySize)

#define Uni_Structured_Grid_Basic_COMPUTE_GLOBAL_INDEX_3D(x, y, z, xSize, ySize) \
  (x) + (y) * (xSize) + (z) * (xSize) * (ySize)

namespace Uni {
namespace StructuredGrid {
namespace Basic {
//

struct DummyForThisGlobalMultiIndex {
  using CellAccessorType = void;
  using GridType         = void;
};

class Index {
public:
  Index(int const& value = 0) : _indent(1),
    _add(1),
    _value(value),
    _global(&_value) {}

  Index(int const& indent,
        int const& add,
        int const& value) : _indent(indent),
    _add(add),
    _value(value),
    _global(&_value) {}

  Index(int const& indent,
        int const& add,
        int const& value,
        int*       global)
    : _indent(indent),
    _add(add),
    _value(value),
    _global(global) {}

  Index(Index const& other) : _indent(other._indent),
    _add(other._add),
    _value(other._value),
    _global(other._global) {}

  ~Index() {}

  void
  initialize(int const& value = 0) {
    _indent = 1;
    _add    = 1;
    _value  = value;
    _global = &_value;
  }

  void
  initialize(int const& indent,
             int const& add,
             int const& value) {
    _indent = indent;
    _add    = add;
    _value  = value;
    _global = &_value;
  }

  void
  initialize(int const& indent,
             int const& add,
             int const& value,
             int*       global) {
    _indent = indent;
    _add    = add;
    _value  = value;
    _global = global;
  }

  int&
  add() {
    return _add;
  }

  int const&
  add() const {
    return _add;
  }

  int&
  value() {
    return _value;
  }

  int const&
  value() const {
    return _value;
  }

  int*
  global() const {
    return _global;
  }

  void
  global(int* global) {
    _global = global;
  }

  Index&
  operator=(Index const& other) {
    _indent = other._indent;
    _add    = other._add;
    _value  = other._value;
    _global = other._global;

    return *this;
  }

  template <typename U>
  Index&
  operator=(U const& value) {
    _value = value;

    return *this;
  }

  operator int() {
    return _value;
  }

  bool
  operator==(Index const& other) const {
    return _value == other._value;
  }

  template <typename U>
  Index&
  operator+=(U const& value) {
    _value   += value;
    *_global += value * _add;

    return *this;
  }

  Index&
  operator++() {
    ++_value;
    *_global += _indent;

    return *this;
  }

  Index
  operator++(int) {
    auto temp = *this;
    ++_value;
    *_global += _indent;

    return temp;
  }

  template <typename U>
  Index&
  operator-=(U const& value) {
    _value   -= value;
    *_global -= value * _add;

    return *this;
  }

  Index&
  operator--() {
    --_value;
    *_global -= _indent;

    return *this;
  }

  Index
  operator--(int) {
    auto temp = *this;
    --_value;
    *_global -= _indent;

    return temp;
  }

private:
  int  _indent;
  int  _add;
  int  _value;
  int* _global;
};

/**
 * \brief
 * Class, which represents a storage (state) of cell indices in each of three
 * dimensions and global (one-dimensional) cell index.
 *
 * \details
 * The class implements MultiIndex interface. It stores a three-dimensional
 *  index and global (one-dimensional) index.
 *
 *
 * \todo
 * 1. Rename to MultiIndexWithGlobalIndex
 *
 * \author Viacheslav Mikerov <SlavaMikerov@gmail.com>
 * \copyright GNU Public License.
 * \version 1.0 \date 2014
 */
template <typename TDerived,
          typename TGrid,
          unsigned TDimensions,
          typename TDerivedTraits = DummyForThisGlobalMultiIndex>
class GlobalMultiIndex
  : public StructuredGrid::MultiIndex<
      typename boost::mpl::if_<
        boost::is_same<TDerived, DummyForThisGlobalMultiIndex>,
        GlobalMultiIndex<TDerived, TGrid, TDimensions, TDerivedTraits>,
        typename boost::mpl::if_<boost::is_same<TDerived, void>,
                                 typename TDerivedTraits::CellAccessorType,
                                 TDerived>::type
        >::type,
      TDimensions> {
public:
  using Derived
          = typename boost::mpl::if_<
          boost::is_same<TDerived, DummyForThisGlobalMultiIndex>,
          GlobalMultiIndex<TDerived, TGrid, TDimensions, TDerivedTraits>,
          typename boost::mpl::if_<boost::is_same<TDerived, void>,
                                   typename TDerivedTraits::CellAccessorType,
                                   TDerived>::type
          >::type;

  using Base = StructuredGrid::MultiIndex<Derived, TDimensions>;

  using GridType = typename boost::mpl::if_<
          boost::is_same<TGrid, void>,
          typename TDerivedTraits::GridType,
          TGrid>::type;

  enum {
    Dimensions = Base::Dimensions
  };

  using VectorDi = typename Base::VectorDi;

public:
  GlobalMultiIndex(GridType const* grid) : _grid(grid) {
    _indices[0].initialize(1, 1, 0);
    int factor = 1;

    for (int d = 1; d < Dimensions; ++d) {
      _indices[d].initialize(_grid->indentSize(d - 1) * factor,
                             _grid->size(d - 1) * factor,
                             0);
      factor *= _grid->size(d - 1);
    }
    _indices[Dimensions].initialize(0);

    for (int i = 0; i < Dimensions; ++i) {
      _indices[i].global(&_indices[Dimensions].value());
    }
  }

  GlobalMultiIndex(GridType const* grid,
                   VectorDi const& index_) : _grid(grid) {
    _indices[0].initialize(1, 1, index_(0));
    int global = index_(0);
    int factor = 1;

    for (int d = 1; d < Dimensions; ++d) {
      _indices[d].initialize(_grid->indentSize(d - 1) * factor,
                             _grid->size(d - 1) * factor,
                             index_(d));
      factor *= _grid->size(d - 1);
      global += factor * index_(d);
    }

    _indices[Dimensions].initialize(global);

    for (int i = 0; i < Dimensions; ++i) {
      _indices[i].global(&_indices[Dimensions].value());
    }
  }

  GlobalMultiIndex(GlobalMultiIndex const& other)
    : _grid(other._grid) {
    _indices[Dimensions].initialize(other._indices[Dimensions].value());

    for (int i = 0; i < Dimensions; ++i) {
      _indices[i] = other._indices[i];
      _indices[i].global(&_indices[Dimensions].value());
    }
  }

  ~GlobalMultiIndex() {}

  GlobalMultiIndex&
  operator=(GlobalMultiIndex const& other) {
    _grid = other._grid;
    _indices[Dimensions].initialize(other._indices[Dimensions].value());

    for (int i = 0; i < Dimensions; ++i) {
      _indices[i] = other._indices[i];
      _indices[i].global(&_indices[Dimensions].value());
    }

    return *this;
  }

  void
  initialize(unsigned const& globalIndex_) {
    std::div_t divt;
    divt.quot = globalIndex_;

    for (unsigned d = 0; d < (Dimensions - 1); ++d) {
      divt        = std::div(divt.quot, _grid->size(d));
      _indices[d] = divt.rem;
    }

    _indices[Dimensions - 1] = divt.quot;

    _indices[Dimensions] = globalIndex_;
  }

  void
  initialize(VectorDi const& index_) {
    int global = 0;

    for (int d = Dimensions - 1; d >= 1; --d) {
      _indices[d] = index_(d);
      global     += index_(d) * _grid->size(d - 1);
    }

    _indices[0] = index_(0);
    global     += index_(0);

    _indices[Dimensions] = global;
  }

  bool
  equal(Derived const& other) const {
    return _indices[Dimensions] == other._indices[Dimensions];
  }

  VectorDi
  index() const {
    return this->Base::index();
  }

  Index&
  index(int const& i) {
    return _indices[i];
  }

  Index const&
  index(int const& i) const {
    return _indices[i];
  }

  int const&
  indexValue(int const& i) const {
    return _indices[i].value();
  }

  int
  globalIndex() const {
    return _indices[Dimensions].value();
  }

  int
  absoluteGlobalIndex(VectorDi const& index_) const {
    int global = 0;

    for (int d = Dimensions - 1; d >= 1; --d) {
      global += index_(d) * _grid->size(d - 1);
    }

    global += index_(0);

    return global;
  }

  int
  relativeGlobalIndex(int const& dimension,
                      int const& offset) const {
    int index = globalIndex();
    index += offset * _indices[dimension].add();

    return index;
  }

  int
  relativeGlobalIndex(int const& dimension,
                      int const& offset,
                      int const& dimension2,
                      int const& offset2) const {
    int index = globalIndex();
    index += offset * _indices[dimension].add();
    index += offset2 * _indices[dimension2].add();

    return index;
  }

  int
  relativeGlobalIndex(VectorDi const& offset) const {
    int index = globalIndex();

    for (int d = 0; d < Dimensions; ++d) {
      index += offset(d) * _indices[d].add();
    }

    return index;
  }

private:
  GridType const*                   _grid;
  std::array<Index, Dimensions + 1> _indices;

  Derived const*
  derived() const {
    return static_cast<Derived const*>(this);
  }

  Derived*
  derived() {
    return static_cast<Derived*>(this);
  }
};
}
}
}

#endif
