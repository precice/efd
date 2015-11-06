#pragma once

#include "BodyForce/functions.hpp"

#include "Simulation/Grid.hpp"
#include "Simulation/IfsfdCellAccessor.hpp"
#include "Simulation/IfsfdMemory.hpp"
#include "Simulation/Private/mpigenerics.hpp"
#include "Simulation/Reporter.hpp"
#include "Simulation/SfsfdCellAccessor.hpp"
#include "Simulation/SfsfdMemory.hpp"

#include "functions.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/ExecutionControl/exception>
#include <Uni/Helpers/macros>

#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

#include <array>
#include <map>

namespace Fluid {
namespace Simulation {
namespace ImmersedBoundary {
template <typename TValue>
class IteratorBackEnd {
  virtual
  ~IteratorBackEnd() {}

  virtual bool
  equals(IteratorBackEnd const* other) = 0;

  virtual TValue&
  dereference() const = 0;

  virtual TValue const&
  dereference() = 0;

  virtual void
  increment() = 0;

  virtual void
  decrement() = 0;
};

template <typename TValue>
class Iterator :
  public Uni::IteratorFacade<Iterator<TValue>, TValue> {
public:
  Iterator(std::unique_ptr<TValue> iterable)
    : _iterable(iterable) {}

  bool
  equals(Iterator const& other) const {
    return _iterable->equals(other._iterable.get());
  }

  TValue const&
  dereference() const {
    return _iterable->dereference();
  }

  void
  increment() {
    _iterable->increment();
  }

  void
  decrement() {
    _iterable->decrement();
  }

private:
  std::unique_ptr < IteratorBackEnd < TValue >> _iterable;
};

template <typename TVector>
class InterfaceCell {
public:
  using Vector = TVector;

  using Scalar = typename Vector::Scalar;

  enum {
    Dimensions = Vector::RowsAtCompileTime
  };

  virtual ~InterfaceCell() {}

  virtual Vector&
  data() = 0;

  virtual Scalar&
  data(unsigned const& dimension) = 0;
};

template <typename TValue>
class IterableBackEnd {
public:
  virtual ~IterableBackEnd() {}

  virtual unsigned
  size() = 0;

  virtual Iterator<TValue> begin() = 0;

  virtual Iterator<TValue> end() = 0;

  virtual Iterator<TValue> find() = 0;
};

template <typename TValue>
class Iterable {
public:
  using BackEnd = std::unique_ptr < IterableBackEnd < TValue >>;

  Iterable(BackEnd const & back_end) : _backEnd(std::move(back_end)) {}

  unsigned
  size() {
    return _backEnd->size();
  }

  Iterator<TValue> begin() {
    return _backEnd->begin();
  }

  Iterator<TValue> end() {
    return _backEnd->begin();
  }

  Iterator<TValue> find() {
    return _backEnd->find();
  }

private:
  BackEnd _backEnd;
};

template <typename TSolverTraits>
class Controller {
public:
  using VectorDsType = typename TSolverTraits::VectorDsType;

  using InterfaceCellType = InterfaceCell<VectorDsType>;

  using IterableType = Iterable<InterfaceCellType>;

  virtual ~Controller() {}

  virtual void
  locateInterfaceCells() = 0;

  virtual void
  precompute() = 0;

  virtual void
  processVelocities() = 0;

  virtual void
  processForces() = 0;

  virtual IterableType
  getVelocityIterable() = 0;

  virtual IterableType
  getForceIterable() = 0;
};
}
}
}
