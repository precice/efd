#ifndef Uni_StructuredGrid_MultiIndex_hpp
#define Uni_StructuredGrid_MultiIndex_hpp

#include <Eigen/Core>

#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>

namespace Uni {
namespace StructuredGrid {
struct DummyForThisMultiIndex {
  using CellAccessorType = void;
  using GridType         = void;
};
/**
 * \brief
 * Abstract template class, which represents a storage (state) of cell indices.
 *
 * \details
 * The class is intended to provide indices for a grid to store a cell position
 * in it. That class is an iterface for a grid iterator which methods are used
 * to implement the iterator functional.
 *
 * \author Viacheslav Mikerov <SlavaMikerov@gmail.com>
 * \copyright GNU Public License.
 * \version 1.0 \date 2014
 */
template <typename TDerived,
          unsigned TDimensions,
          typename TDerivedTraits = DummyForThisMultiIndex>
class MultiIndex {
public:
  enum {
    Dimensions = TDimensions
  };

  using Derived = typename boost::mpl::if_<
          boost::is_same<TDerived, DummyForThisMultiIndex>,
          MultiIndex,
          typename boost::mpl::if_<boost::is_same<TDerived, void>,
                                   typename TDerivedTraits::CellAccessorType,
                                   TDerived>::type
          >::type;

  using VectorDi = Eigen::Matrix<int, Dimensions, 1>;

public:
  /**
   * The followin methods must be implemented
   */
  /*
   *  virtual void
   *  initialize(int const&,
   *          int const&,
   *          int const&) = 0;
   *
   *  virtual bool
   *  equal(MultiIndex const&) const {
   *  return
   *  }
   *
   *  virtual T&
   *  index(int const&) = 0;
   *
   *  virtual T const&
   *  index(int const&) const = 0;
   *
   *  virtual int const&
   *  indexValue(int const&) const = 0;
   */
  VectorDi
  index() const {
    return indexValues();
  }

  VectorDi
  indexValues() const {
    VectorDi result;

    for (int i = 0; i < Dimensions; ++i) {
      result(i) = derived().indexValue(i);
    }

    return result;
  }

  void
  operator()(VectorDi const& i) {
    derived().initialize(i);
  }

  bool
  operator==(Derived const& other) const {
    return derived().equal(other);
  }

  int const&
  operator()(int const& i) const { return derived().indexValue(i); }

  int&
  operator[](int const& i) { return derived().indexValue(i); }

  int const&
  operator[](int const& i) const { return derived().indexValue(i); }

  VectorDi
  operator()() const {
    return derived().indexValues();
  }

  void
  move(int const& dimension,
       int const& offset) {
    derived().index(dimension) += offset;
  }

  void
  move(int const& dimension,
       int const& offset,
       int const& dimension2,
       int const& offset2) {
    derived().index(dimension)  += offset;
    derived().index(dimension2) += offset2;
  }

  void
  move(VectorDi const& i) {
    for (int d = 0; d < Dimensions; ++d) {
      derived().index(d) += i(d);
    }
  }

  Derived
  absolute(VectorDi const& i) const {
    Derived result(derived());
    result.initialize(i);

    return result;
  }

  Derived
  relative(int const& dimension,
           int const& offset) const {
    Derived result(derived());
    result.move(dimension, offset);

    return result;
  }

  Derived
  relative(int const& dimension,
           int const& offset,
           int const& dimension2,
           int const& offset2) const {
    Derived result(derived());
    result.move(dimension, offset,
                dimension2, offset2);

    return result;
  }

  Derived
  relative(VectorDi const& i) const {
    Derived result(derived());
    result.move(i);

    return result;
  }

  VectorDi
  absoluteIndex(VectorDi const& i) const {
    return i;
  }

  VectorDi
  relativeIndex(int const& dimension,
                int const& offset) const {
    VectorDi result(indexValues());
    result(dimension) += offset;

    return result;
  }

  VectorDi
  relativeIndex(int const& dimension,
                int const& offset,
                int const& dimension2,
                int const& offset2) const {
    VectorDi result(indexValues());
    result(dimension)  += offset;
    result(dimension2) += offset2;

    return result;
  }

  VectorDi
  relativeIndex(VectorDi const& i) const {
    return (indexValues() + i).eval();
  }

private:
  Derived const&
  derived() const {
    return *static_cast<Derived const*>(this);
  }

  Derived&
  derived() {
    return *static_cast<Derived*>(this);
  }
};
}
}
#endif
