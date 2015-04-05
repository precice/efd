#pragma once

#include <precice/SolverInterface.hpp>

#include <Uni/Helpers/macros>

#include <Eigen/Core>

namespace Structure {
class Simulation {
public:
  using VectorDs = Eigen::Matrix<double, 2, 1>;

public:
  Simulation() {}

  Simulation(Simulation const&) = delete;

  ~Simulation() {}

  Simulation&
  operator=(Simulation const&) = delete;

  void
  initialize(precice::SolverInterface* precice_interface) {
    // _type = type;
    _preciceInterface = precice_interface;
  }

  bool
  iterate() {
    return false;
    //
  }

  Uni_PublicProperty(int,      type)
  Uni_PublicProperty(VectorDs, velocity)

private:
  precice::SolverInterface* _preciceInterface;
};
}
