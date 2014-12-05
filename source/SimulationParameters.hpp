#ifndef  FsiSimulation_SimulationParameters_hpp
#define  FsiSimulation_SimulationParameters_hpp

#include <Eigen/Core>

namespace FsiSimulation {
template <typename Scalar, int D>
class SimulationParameters {
public:
  typedef Eigen::Matrix<Scalar, D, 1> VectorDs;

public:
  SimulationParameters() {}

  SimulationParameters(SimulationParameters const& other) = delete;

  ~SimulationParameters() {}

  SimulationParameters&
  operator=(SimulationParameters const& other) = delete;

  Scalar&
  re() { return _re; }

  Scalar const&
  re() const { return _re; }

  Scalar&
  gamma() { return _gamma; }

  Scalar const&
  gamma() const { return _gamma; }

  Scalar&
  tau() { return _tau; }

  Scalar const&
  tau() const { return _tau; }

  VectorDs&
  g() { return _g; }

  VectorDs const&
  g() const { return _g; }

  Scalar&
  g(int const& i) { return _g(i); }

  Scalar const&
  g(int const& i) const { return _g(i); }

private:
  Scalar   _re;
  Scalar   _gamma;
  Scalar   _tau;
  VectorDs _g;
};
}

#endif
