#ifndef FsiSimulation_FluidSimulation_Parameters_hpp
#define FsiSimulation_FluidSimulation_Parameters_hpp

#include <Eigen/Core>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TScalar, int TD>
class Parameters {
public:
  typedef Eigen::Matrix<TScalar, TD, 1> VectorDs;

public:
  Parameters() {}

  Parameters(Parameters const& other) = delete;

  ~Parameters() {}

  Parameters&
  operator=(Parameters const& other) = delete;

  TScalar&
  re() {
    return _re;
  }

  TScalar const&
  re() const {
    return _re;
  }

  TScalar&
  gamma() {
    return _gamma;
  }

  TScalar const&
  gamma() const {
    return _gamma;
  }

  TScalar&
  tau() {
    return _tau;
  }

  TScalar const&
  tau() const {
    return _tau;
  }

  VectorDs&
  g() {
    return _g;
  }

  VectorDs const&
  g() const {
    return _g;
  }

  TScalar&
  g(int const& i) {
    return _g(i);
  }

  TScalar const&
  g(int const& i) const {
    return _g(i);
  }

private:
  TScalar  _re;
  TScalar  _gamma;
  TScalar  _tau;
  VectorDs _g;
};
}
}

#endif
