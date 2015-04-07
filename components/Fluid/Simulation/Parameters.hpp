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
  Parameters() {
    _g = VectorDs::Zero();
  }

  Parameters(Parameters const& other) = delete;

  ~Parameters() {}

  Parameters&
  operator=(Parameters const& other) = delete;

  TScalar const&
  re() const {
    return _re;
  }

  TScalar&
  re() {
    return _re;
  }

  TScalar const&
  gamma() const {
    return _gamma;
  }

  TScalar&
  gamma() {
    return _gamma;
  }

  TScalar const&
  tau() const {
    return _tau;
  }

  TScalar&
  tau() {
    return _tau;
  }

  TScalar const&
  alpha() const {
    return _alpha;
  }

  TScalar&
  alpha() {
    return _alpha;
  }

  VectorDs const&
  g() const {
    return _g;
  }

  VectorDs&
  g() {
    return _g;
  }

  TScalar const&
  g(int const& i) const {
    return _g(i);
  }

  TScalar&
  g(int const& i) {
    return _g(i);
  }

private:
  TScalar  _re;
  TScalar  _gamma;
  TScalar  _tau;
  TScalar  _alpha;
  VectorDs _g;
};
}
}

#endif
