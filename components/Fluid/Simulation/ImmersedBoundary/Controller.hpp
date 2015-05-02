#pragma once

#include <Uni/IteratorFacade>

#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
template <typename TValue>
class IteratorBackEnd {
public:
  virtual
  ~IteratorBackEnd() {}

  virtual std::unique_ptr<IteratorBackEnd>
  clone() const = 0;

  virtual bool
  equals(std::unique_ptr<IteratorBackEnd> const& other) = 0;

  virtual TValue const&
  dereference() const = 0;

  virtual TValue&
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
  using BackEndType = IteratorBackEnd<TValue>;

  using UniqueBackEndType = std::unique_ptr<BackEndType>;

  Iterator(UniqueBackEndType&& back_end)
    : _backEnd(std::move(back_end)) {}

  Iterator(Iterator const& other)
    : _backEnd(std::move(other._backEnd->clone())) {}

  Iterator&
  operator=(Iterator const& other) {
    _backEnd = std::move(other._backEnd->clone());

    return *this;
  }

  bool
  equals(Iterator const& other) const {
    return _backEnd->equals(other._backEnd);
  }

  TValue&
  dereference() const {
    return const_cast<TValue&>(_backEnd->dereference());
  }

  void
  increment() {
    _backEnd->increment();
  }

  void
  decrement() {
    _backEnd->decrement();
  }

private:
  UniqueBackEndType _backEnd;
};

template <typename TValue>
class IterableBackEnd {
public:
  using Value = TValue;

  using IteratorType = Iterator<Value>;

  virtual ~IterableBackEnd() {}

  virtual unsigned
  size() = 0;

  virtual IteratorType
  begin() = 0;

  virtual IteratorType
  end() = 0;

  virtual IteratorType
  find(unsigned const& global_index) = 0;
};

template <typename TValue>
class CompoundIteratorBackEnd : public IteratorBackEnd<TValue> {
public:
  using Base = IteratorBackEnd<TValue>;

  using IteratorType = Iterator<TValue>;

  using IterableBackEndType = IterableBackEnd<TValue>;

  CompoundIteratorBackEnd(
    IteratorType const&                         it,
    bool const&                                 is_first,
    std::shared_ptr<IterableBackEndType> const& first_iterable,
    std::shared_ptr<IterableBackEndType> const& second_iterable)
    : _it(it),
    _isFirst(is_first),
    _firstIterable(first_iterable),
    _secondIterable(second_iterable) {}

  CompoundIteratorBackEnd(CompoundIteratorBackEnd&) = delete;

  virtual
  ~CompoundIteratorBackEnd() {}

  CompoundIteratorBackEnd&
  operator=(CompoundIteratorBackEnd&) = delete;

  virtual std::unique_ptr<Base>
  clone() const {
    return std::unique_ptr<Base>(
      new CompoundIteratorBackEnd(_it,
                                  _isFirst,
                                  _firstIterable,
                                  _secondIterable));
  }

  virtual bool
  equals(std::unique_ptr<Base> const& other) {
    return _it == reinterpret_cast<CompoundIteratorBackEnd*>
           (other.get())->_it
           && _isFirst == reinterpret_cast<CompoundIteratorBackEnd*>
           (other.get())->_isFirst;
  }

  virtual TValue const&
  dereference() const {
    return *_it;
  }

  virtual TValue&
  dereference() {
    return *_it;
  }

  virtual void
  increment() {
    if (_it == _secondIterable->end()) {
      return;
    }

    ++_it;

    if (_it == _firstIterable->end()) {
      _isFirst = false;
      _it = _secondIterable->begin();
    }
  }

  virtual void
  decrement() {
    if (_it == _firstIterable->begin()) {
      return;
    }

    if (_it == _secondIterable->begin()) {
      _isFirst = true;
      _it = _firstIterable->end();
    }
    --_it;
  }

private:
  IteratorType _it;
  bool         _isFirst;

  std::shared_ptr<IterableBackEndType> _firstIterable;
  std::shared_ptr<IterableBackEndType> _secondIterable;
};

template <typename TValue>
class CompoundIterableBackEnd : public IterableBackEnd<TValue> {
public:
  using Value = TValue;

  using Base = IterableBackEnd<Value>;

  using IteratorType = typename Base::IteratorType;

  using IteratorBackEndType = typename IteratorType::BackEndType;

  using UniqueIteratorBackEndType = typename IteratorType::UniqueBackEndType;

  using IterableBackEndType = Base;

  using SharedIterableBackEndType = std::shared_ptr<IterableBackEndType>;

  CompoundIterableBackEnd(
    SharedIterableBackEndType const& first_iterable,
    SharedIterableBackEndType const& second_iterable)
    : _firstIterable(first_iterable),
    _secondIterable(second_iterable) {}

  virtual ~CompoundIterableBackEnd() {}

  virtual unsigned
  size() {
    return _firstIterable->size() + _secondIterable->size();
  }

  virtual IteratorType
  begin() {
    if (_firstIterable->begin() == _firstIterable->end()) {
      return IteratorType(
        UniqueIteratorBackEndType(
          new CompoundIteratorBackEnd<Value>(
            _secondIterable->begin(),
            false,
            _firstIterable,
            _secondIterable)));
    }

    return IteratorType(
      UniqueIteratorBackEndType(
        new CompoundIteratorBackEnd<Value>(
          _firstIterable->begin(),
          true,
          _firstIterable,
          _secondIterable)));
  }

  virtual IteratorType
  end() {
    return IteratorType(
      UniqueIteratorBackEndType(
        new CompoundIteratorBackEnd<Value>(
          _secondIterable->end(),
          false,
          _firstIterable,
          _secondIterable)));
  }

  virtual IteratorType
  find(unsigned const& global_index) {
    auto find_it = _firstIterable->find(global_index);

    if (find_it == _firstIterable->end()) {
      find_it = _secondIterable->find(global_index);

      return IteratorType(
        UniqueIteratorBackEndType(
          new CompoundIteratorBackEnd<Value>(
            find_it,
            false,
            _firstIterable,
            _secondIterable)));
    } else {
      return IteratorType(
        UniqueIteratorBackEndType(
          new CompoundIteratorBackEnd<Value>(
            find_it,
            true,
            _firstIterable,
            _secondIterable)));
    }
  }

private:
  SharedIterableBackEndType _firstIterable;
  SharedIterableBackEndType _secondIterable;
};

template <typename TValue>
class Iterable {
public:
  using BackEndType = IterableBackEnd<TValue>;

  using SharedBackEndType = std::shared_ptr<BackEndType>;

  Iterable(SharedBackEndType const& back_end)
    : _backEnd(back_end) {}

  Iterable&
  operator=(Iterable const& other) {
    _backEnd = other._backEnd;

    return *this;
  }

  unsigned
  size() {
    return _backEnd->size();
  }

  Iterator<TValue> begin() {
    return _backEnd->begin();
  }

  Iterator<TValue> end() {
    return _backEnd->end();
  }

  Iterator<TValue> find(unsigned const& global_index) {
    return _backEnd->find(global_index);
  }

private:
  SharedBackEndType _backEnd;
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

  InterfaceCell&
  operator=(InterfaceCell&) = delete;

  virtual Vector const&
  data() const = 0;

  virtual Vector&
  data() = 0;

  virtual Scalar const&
  data(unsigned const& dimension) const = 0;

  virtual Scalar&
  data(unsigned const& dimension) = 0;

  virtual unsigned
  globalIndex() const = 0;
};

template <typename TVector>
class Controller {
public:
  using VectorType = TVector;

  using InterfaceCellType = InterfaceCell<VectorType>;

  using IterableType = Iterable<InterfaceCellType>;

  virtual ~Controller() {}

  virtual void
  initialize() = 0;

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
