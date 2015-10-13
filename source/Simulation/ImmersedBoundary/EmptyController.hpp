#pragma once

#include "Controller.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
template <typename TVector>
class EmptyInterfaceCell : public InterfaceCell<TVector> {
public:
  using Vector = TVector;

  using Scalar = typename Vector::Scalar;

  ~EmptyInterfaceCell() {}

  Vector const&
  data() const {
    static Vector dummy;

    return dummy;
  }

  Vector&
  data() {
    static Vector dummy;

    return dummy;
  }

  Scalar const&
  data(unsigned const& dimension) const {
    ((void)dimension);
    static Scalar dummy;

    return dummy;
  }

  Scalar&
  data(unsigned const& dimension) {
    ((void)dimension);
    static Scalar dummy;

    return dummy;
  }

  unsigned
  globalIndex() const {
    return 0;
  }
};

template <typename TVector>
class EmptyIteratorBackEnd : public IteratorBackEnd < InterfaceCell < TVector >>
{
public:
  using Value = InterfaceCell<TVector>;

  using Base = IteratorBackEnd<Value>;

  ~EmptyIteratorBackEnd() {}

  std::unique_ptr<Base>
  clone() const {
    return std::unique_ptr<Base>(new EmptyIteratorBackEnd());
  }

  bool
  equals(std::unique_ptr<Base> const& other) {
    ((void)other);
    return true;
  }

  Value&
  dereference() const {
    static EmptyInterfaceCell<TVector> dummy;

    return dummy;
  }

  void
  increment() {}

  void
  decrement() {}
};

template <typename TVector>
class EmptyIterableBackEnd :
  public IterableBackEnd < InterfaceCell < TVector >> {
public:
  using Vector = TVector;

  using Value = InterfaceCell<Vector>;

  using Base = IterableBackEnd<Value>;

  using BaseIteratorType = typename Base::IteratorType;

  using BaseUniqueIteratorBackEndType
          = typename BaseIteratorType::UniqueBackEndType;

  using IteratorBackEndType = EmptyIteratorBackEnd<TVector>;

  ~EmptyIterableBackEnd() {}

  unsigned
  size() {
    return 0u;
  }

  BaseIteratorType
  begin() {
    return BaseIteratorType(
      BaseUniqueIteratorBackEndType(new IteratorBackEndType));
  }

  BaseIteratorType
  end() {
    return BaseIteratorType(
      BaseUniqueIteratorBackEndType(new IteratorBackEndType));
  }

  BaseIteratorType
  find(unsigned const& global_index) {
    ((void)global_index);
    return BaseIteratorType(
      BaseUniqueIteratorBackEndType(new IteratorBackEndType));
  }
};

template <typename TVector>
class EmptyController : public Controller<TVector> {
public:
  using VectorType = TVector;

  using BaseType = Controller<VectorType>;

  using IterableType = typename BaseType::IterableType;

private:
  using BaseInterfaceCellType = typename BaseType::InterfaceCellType;

  using InterfaceCellType = EmptyInterfaceCell<VectorType>;

  using BaseIterableBackEndType
          = typename IterableType::BackEndType;

  using BaseSharedIterableBackEndType
          = typename IterableType::SharedBackEndType;

  using IterableBackEndType = EmptyIterableBackEnd<VectorType>;

public:
  ~EmptyController() {}

  void
  initialize(precice::SolverInterface* precice_interface) {
    ((void)precice_interface);
  }

  void
  precompute() {}

  void
  processVelocities() {}

  void
  processForces() {}

  IterableType
  getVelocityIterable() {
    return IterableType(
      BaseSharedIterableBackEndType(new IterableBackEndType));
  }

  IterableType
  getForceIterable() {
    return IterableType(
      BaseSharedIterableBackEndType(new IterableBackEndType));
  }
};
}
}
}
