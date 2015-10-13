#ifndef Uni_StructuredGrid_Basic_Grid_hpp
#define Uni_StructuredGrid_Basic_Grid_hpp

#include <Uni/StructuredGrid/GridIterator>

#include <Eigen/Core>

#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>

#include <array>

namespace Uni {
namespace StructuredGrid {
namespace Basic {
struct DummyForThisGrid {
  using CellAccessorType = void;
  using GridType         = void;
};
/**
 * \brief
 * Abstract class, which represents an entity of a three-dimensional structured
 * grid. Implements some methods of Grid interface, which most probably
 * would be considered as derivative methods of other methods.
 *
 * \details
 * The class is intended to be a comfortable form of the Grid interface in
 * some situations since it implements some methods of the iterface based on
 * other methods of the interface.
 *
 * \author Viacheslav Mikerov <SlavaMikerov@gmail.com>
 * \copyright GNU Public License.
 * \version 1.0 \date 2014
 */
template <typename TDerived,
          typename TCellAccessor,
          unsigned TDimensions,
          typename TDerivedTraits = DummyForThisGrid>
class Grid {
public:
  using Derived
          = typename boost::mpl::if_<
          boost::is_same<TDerived, DummyForThisGrid>,
          Grid,
          typename boost::mpl::if_<boost::is_same<TDerived, void>,
                                   typename TDerivedTraits::GridType,
                                   TDerived>::type
          >::type;

  using CellAccessor = typename boost::mpl::if_<
          boost::is_same<TCellAccessor, void>,
          typename TDerivedTraits::CellAccessorType,
          TCellAccessor>::type;

  enum {
    Dimensions = TDimensions
  };

  using Iterator
          = StructuredGrid::GridIterator<Derived, CellAccessor, TDimensions>;

  using VectorDi = Eigen::Matrix<int, Dimensions, 1>;

  using Factory = typename Iterator::Factory;

public:
  Grid() {}

  Grid(Grid const& other) : _factory(other._factory),
    _size(other._size),
    _indent(other._indent),
    _innerLimit(other._innerLimit),
    _innerSize(other._innerSize),
    _indentSize(other._indentSize) {}

  ~Grid() {}

  Grid&
  operator=(Grid const& other) {
    _factory    = other._factory;
    _size       = other._size;
    _indent     = other._indent;
    _innerLimit = other._innerLimit;
    _innerSize  = other._innerSize;
    _indentSize = other._indentSize;

    return *this;
  }

  void
  initialize(VectorDi const& size,
             VectorDi const& leftIndent = VectorDi::Zero(),
             VectorDi const& rightIndent = VectorDi::Zero(),
             Factory const&  factory = defaultFactory()) {
    _factory    = factory;
    _size       = size;
    _indent[0]  = leftIndent;
    _indent[1]  = rightIndent;
    _innerLimit = _size - _indent[1];
    _innerSize  = _innerLimit - _indent[0];
    _indentSize = _indent[0] + _indent[1];
  }

  void
  setIndents(VectorDi const& leftIndent = VectorDi::Zeros(),
             VectorDi const& rightIndent = VectorDi::Zeros()) {
    _indent[0]  = leftIndent;
    _indent[1]  = rightIndent;
    _innerLimit = _size - _indent[1];
    _innerSize  = _innerLimit - _indent[0];
    _indentSize = _indent[0] + _indent[1];
  }

  Iterator
  begin() const {
    return Iterator(derived(), _factory);
  }

  Iterator
  end() const {
    return Iterator(derived(), true, _factory);
  }

  Iterator
  at(VectorDi const& index) const {
    return Iterator(derived(), index, _factory);
  }

  VectorDi const&
  size() const {
    return _size;
  }

  int const&
  size(int const& dimension) const {
    return _size(dimension);
  }

  VectorDi const&
  indent(int const& i) const {
    return _indent[i];
  }

  VectorDi const&
  leftIndent() const {
    return _indent[0];
  }

  int const&
  leftIndent(int const& dimension) const {
    return _indent[0](dimension);
  }

  VectorDi const&
  rightIndent() const {
    return _indent[1];
  }

  int const&
  rightIndent(int const& dimension) const {
    return _indent[1](dimension);
  }

  VectorDi const&
  innerLimit() const {
    return _innerLimit;
  }

  int const&
  innerLimit(int const& dimension) const {
    return _innerLimit(dimension);
  }

  VectorDi const&
  innerSize() const {
    return _innerSize;
  }

  int const&
  innerSize(int const& dimension) const {
    return _innerSize(dimension);
  }

  VectorDi const&
  indentSize() const {
    return _indentSize;
  }

  int const&
  indentSize(int const& dimension) const {
    return _indentSize(dimension);
  }

  static Factory
  defaultFactory() {
    static Factory defaultFactory = [] (Grid const* grid, VectorDi const& i) {
                                      return CellAccessor(grid, i);
                                    };

    return defaultFactory;
  }

protected:
  Factory                 _factory;
  VectorDi                _size;
  std::array<VectorDi, 2> _indent;
  VectorDi                _innerLimit;
  VectorDi                _innerSize;
  VectorDi                _indentSize;

private:
  Derived const*
  derived() const {
    return reinterpret_cast<Derived const*>(this);
  }

  Derived*
  derived() {
    return reinterpret_cast<Derived*>(this);
  }
};
}
}
}
#endif
