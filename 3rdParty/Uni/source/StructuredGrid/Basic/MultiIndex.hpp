#ifndef Uni_StructuredGrid_Basic_MultiIndex_hpp
#define Uni_StructuredGrid_Basic_MultiIndex_hpp

#include <Uni/StructuredGrid/MultiIndex>

namespace Uni {
namespace StructuredGrid {
namespace Basic {
template <typename TDerivedTraits>
class MultiIndex : public StructuredGrid::MultiIndex<TDerivedTraits> {
public:
  using Base = StructuredGrid::MultiIndex<TDerivedTraits>;

  enum {
    Dimensions = Base::Dimensions
  };

  typedef typename Base::VectorDi VectorDi;

public:
  MultiIndex() {}

  MultiIndex(VectorDi const& i) : _index(i) {}

  MultiIndex(MultiIndex const& other) : _index(other._index) {}

  ~MultiIndex() {}

  MultiIndex&
  operator=(MultiIndex const& other) {
    _index = other._index;

    return *this;
  }

  VectorDi
  indexValues() const {
    return _index;
  }

  void
  initialize(VectorDi const& i) {
    _index = i;
  }

  bool
  equal(MultiIndex const& other) const {
    return _index == other._index;
  }

  int&
  index(int const& i) {
    return *(_index.data() + i);
  }

  int const&
  index(int const& i) const {
    return *(_index.data() + i);
  }

  int const&
  indexValue(int const& i) const {
    return *(_index.data() + i);
  }

private:
  VectorDi _index;
};
}
}
}

#endif
