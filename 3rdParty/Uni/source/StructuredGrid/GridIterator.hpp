#ifndef Uni_StructuredGrid_GridIterator_hpp
#define Uni_StructuredGrid_GridIterator_hpp

#include <Eigen/Core>

#include <boost/iterator/iterator_facade.hpp>

#include <functional>

namespace Uni {
namespace StructuredGrid {
/**
 * \brief
 * Template class, which represents a random access iterator (in usual C++
 *  sense) over a three-dimensional grid.
 *
 * \details
 * It inherits the functional part like overloaded operators from IteratorFacade
 * class.
 * In this class only implemented the essential methods for iterator, that then
 * are used by IteratorFacade's methods to expose the functionality.
 *
 * With a template parameter you can decide what is the result of dereference
 * operation by inheritence from CellAccessor class (or your custom class with
 * the same methods).
 *
 * @tparam CellAccessor (defult CellAccessor) is a storage for a
 * three-dimensional
 * index,
 *  which is returned as a value of iterator.
 *
 * \todo
 * 1. Remove global index form the factor.
 *
 * \author Viacheslav Mikerov <SlavaMikerov@gmail.com>
 * \copyright GNU Public License.
 * \version 1.0 \date 2014
 */

template <typename TGrid, typename TCellAccessor, unsigned TDimensions>
class GridIterator
  : public boost::iterator_facade
    <GridIterator<TGrid, TCellAccessor, TDimensions>,
     TCellAccessor,
     std::random_access_iterator_tag,
     TCellAccessor const> {
  friend class boost::iterator_core_access;

public
  :
  using Base = boost::iterator_facade
               <GridIterator<TGrid, TCellAccessor, TDimensions>,
                TCellAccessor,
                std::random_access_iterator_tag,
                TCellAccessor const>;

  using GridType = TGrid;

  using CellAccessor = TCellAccessor;

  enum {
    Dimensions = TDimensions
  };

  typedef Eigen::Matrix<int, Dimensions, 1> VectorDi;

  typedef std::function<CellAccessor(GridType const*, VectorDi const&)> Factory;

  using Difference = typename Base::difference_type;
  using Pointer    = typename Base::pointer;
  using Reference  = typename Base::reference;
  using Value      = typename Base::value_type;

public:
  GridIterator(GridType const* grid,
               Factory const&  factory = defaultFactory())
    : _grid(grid),
    _multiIndex(factory(grid, getBeginningVector()))
  {}

  GridIterator(GridType const* grid,
               VectorDi const& index,
               Factory const&  factory = defaultFactory())
    : _grid(grid),
    _multiIndex(factory(grid, index))
  {}

  GridIterator(GridType const* grid,
               bool const&     isEnd,
               Factory const&  factory = defaultFactory())
    : _grid(grid),
    _multiIndex(factory(grid, getEndVector()))
  { ((void)isEnd); }

  GridIterator(GridIterator const& other)
    : _grid(other._grid),
    _multiIndex(other._multiIndex)
  {}

  GridIterator const&
  operator=(GridIterator const& other) {
    _grid       = other._grid;
    _multiIndex = other._multiIndex;

    return *this;
  }

  static Factory
  defaultFactory() {
    static Factory defaultFactory
      = [] (GridType const* grid, VectorDi const& i) {
          return CellAccessor(grid, i);
        };

    return defaultFactory;
  }

private:
  bool
  equal(GridIterator const& other) const {
    return _multiIndex == other._multiIndex;
  }

  Reference
  dereference() const {
    return _multiIndex;
  }

  void
  increment() {
    if (_multiIndex(Dimensions - 1) < _grid->innerLimit(Dimensions - 1)) {
      ++_multiIndex.index(0);
    } else {
      return;
    }

    for (int i = 0; i < (Dimensions - 1); ++i) {
      if (_multiIndex(i) >= _grid->innerLimit(i)) {
        ++_multiIndex.index(i + 1);
        _multiIndex.index(i) = _grid->leftIndent(i);
      } else {
        break;
      }
    }
  }

  void
  decrement() {
    if (_multiIndex(Dimensions - 1) >= _grid->leftIndent(Dimensions - 1)) {
      --_multiIndex.index(0);
    } else {
      return;
    }

    for (int i = 0; i < (Dimensions - 1); ++i) {
      if (_multiIndex(i) < _grid->leftIndent(i)) {
        --_multiIndex.index(i + 1);
        _multiIndex.index(i) = _grid->innerLimit(i) - 1;
      } else {
        break;
      }
    }
  }

  void
  advance(int const& dist) {
    bool isPositive   = true;
    int  unsignedDist = dist;

    if (dist < 0) {
      isPositive   = false;
      unsignedDist = -dist;
    }

    VectorDi distanceVector;

    for (int i = 0; i < (Dimensions - 1); ++i) {
      std::div_t temp = div(unsignedDist, _grid->innerSize(i));
      unsignedDist      = temp.quot;
      distanceVector(i) = temp.rem;
    }
    distanceVector(Dimensions - 1) = unsignedDist;

    if (isPositive) {
      for (int i = 0; i < (Dimensions - 1); ++i) {
        _multiIndex.index(i) += distanceVector(i);

        for (int j = 0; j <= i; ++j) {
          if (_multiIndex(i) >= _grid->innerLimit(i)) {
            ++_multiIndex.index(i + 1);
            _multiIndex.index(i) = _multiIndex(i) - _grid->innerSize(i);
          } else {
            break;
          }
        }
      }

      _multiIndex.index(Dimensions - 1) += distanceVector(Dimensions - 1);

      if (_multiIndex(Dimensions - 1) >= _grid->innerLimit(Dimensions - 1)) {
        VectorDi newValue;

        _multiIndex.initialize(getEndVector());
      }
    } else {
      for (int i = 0; i < (Dimensions - 1); ++i) {
        _multiIndex.index(i) -= distanceVector(i);

        for (int j = 0; j <= i; ++j) {
          if (_multiIndex(i) < _grid->leftIndent(i)) {
            --_multiIndex.index(i + 1);
            _multiIndex.index(i) = _multiIndex(i) +
                                   _grid->innerSize(i) -
                                   _grid->leftIndent(i);
          } else {
            break;
          }
        }
      }

      _multiIndex.index(Dimensions - 1) -= distanceVector(Dimensions - 1);

      if (_multiIndex(Dimensions - 1) < _grid->leftIndent(Dimensions - 1)) {
        _multiIndex.initialize(getBeginningVector());
      }
    }
  }

  int
  distanceTo(GridIterator const& other) const {
    VectorDi distanceVector;

    for (int i = 0; i <= (Dimensions - 1); ++i) {
      distanceVector(i) = other._multiIndex(i) - _multiIndex(i);
    }

    int factor = 1;
    int dist   = 0;

    for (int i = 0; i <= (Dimensions - 1); ++i) {
      dist   += distanceVector(i) * factor;
      factor *= _grid->innerSize(i);
    }

    return dist;
  }

  int
  distance() const {
    VectorDi distanceVector;

    for (int i = 0; i <= (Dimensions - 1); ++i) {
      distanceVector(i) = _multiIndex(i) - _grid->leftIndent(i);
    }

    int factor = 1;
    int dist   = 0;

    for (int i = 0; i <= (Dimensions - 1); ++i) {
      dist   += distanceVector(i) * factor;
      factor *= _grid->innerSize(i);
    }

    return dist;
  }

  VectorDi
  getBeginningVector() {
    return _grid->leftIndent();
  }

  VectorDi
  getEndVector() {
    VectorDi end(_grid->leftIndent());
    end(Dimensions - 1) = _grid->innerLimit(Dimensions - 1);

    return end;
  }

  GridType const* _grid;
  CellAccessor    _multiIndex;
};
}
}

#endif
