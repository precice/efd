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

  int const&
  outerLayerSize() const {
    return _outerLayerSize;
  }

  int&
  outerLayerSize() {
    return _outerLayerSize;
  }

  int const&
  innerLayerSize() const{
    return _innerLayerSize;
  }

  int&
  innerLayerSize() {
    return _innerLayerSize;
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
  int      _outerLayerSize;
  int      _innerLayerSize;
  VectorDs _g;
};
}
}

#endif
