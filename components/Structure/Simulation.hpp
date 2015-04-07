#pragma once

#include <precice/SolverInterface.hpp>

#include <Uni/Helpers/macros>

#include <Eigen/Core>

#include <memory>

namespace Structure {
class Simulation {
public:
  using VectorDd = Eigen::Matrix<double, 2, 1>;
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
//    double dt = _preciceInterface->initialize();

    _meshId   = _preciceInterface->getMeshID("BodyMesh");
    _mesh.reset(new precice::MeshHandle(_preciceInterface->getMeshHandle("BodyMesh")));
    _forcesID = _preciceInterface->getDataID("Forces", _meshId);
  }

  bool
  iterate() {
    if (!_preciceInterface->isCouplingOngoing()) {
      return false;
    }

    static VectorDd velocity = VectorDd::Zero();

    auto it = _mesh->vertices().begin();

    for (std::size_t i = 0; i < _mesh->vertices().size(); ++i, it++) {
      VectorDd coords(it.vertexCoords());
      coords = _shift.cast<double>() + coords;
      _preciceInterface->writeVectorData(_forcesID,
                                         it.vertexID(),
                                         velocity.data());
    }

//    double dt = _preciceInterface->advance(0.0);

    return true;
  }

  Uni_PublicProperty(int,      type)
  Uni_PublicProperty(VectorDs, velocity)

private:
  precice::SolverInterface*            _preciceInterface;
  std::unique_ptr<precice::MeshHandle> _mesh;
  VectorDs                             _shift;
  int                                  _meshId;
  int                                  _forcesID;
};
}
