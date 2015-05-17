#pragma once

#include "Configuration.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/Helpers/macros>

#include <Eigen/Core>

#include <memory>
#include <vector>

namespace Structure {
class Simulation {
public:
  using VectorDd = Eigen::Matrix<long double, 3, 1>;
  using VectorDs = Eigen::Matrix<long double, 3, 1>;

public:
  Simulation(Configuration const* configuration) {
    _dimensions       = configuration->get<unsigned>("/Dimensions");
    _environmentForce = configuration->get < Eigen::Matrix < long double, 3, 1 >>
    ("/EnvironmentForce");
  }

  Simulation(Simulation const&) = delete;

  ~Simulation() {}

  Simulation&
  operator=(Simulation const&) = delete;

  void
  initialize(precice::SolverInterface* precice_interface) {
    // _type = type;
    _preciceInterface = precice_interface;
    _dt               = _preciceInterface->initialize();

    if (!_preciceInterface->hasMesh("BodyMesh")) {
      throwException("Precice does not have 'BodyMesh' in its configuration");
    }

    _meshId = _preciceInterface->getMeshID("BodyMesh");
    // _forcesID = _preciceInterface->getDataID("Forces", _meshId);
    _displacementsID = _preciceInterface->getDataID("Displacements", _meshId);

    unsigned vertices_size = _preciceInterface->getMeshVertexSize(_meshId);

    _vertexIds.resize(vertices_size);
    _displacements.resize(_dimensions * vertices_size);

    for (unsigned i = 0; i < vertices_size; ++i) {
      _vertexIds[i] = i;
    }
    _lastPosition    = VectorDs::Zero();
    _currentPosition = VectorDs::Zero();
  }

  bool
  iterate() {
    if (!_preciceInterface->isCouplingOngoing()) {
      return false;
    }

    // if (_preciceInterface->isActionRequired(
    // precice::constants::actionWriteIterationCheckpoint())) {
    // _preciceInterface->fulfilledAction(
    // precice::constants::actionWriteIterationCheckpoint());
    // }

    // _preciceInterface->resetMesh(_meshId);
    // _preciceInterface->setMeshVertices(_meshId, size, *positions, *ids);
    // _preciceInterface->setMeshEdge(_meshId, firstId, secondId);

    // static VectorDd velocity = VectorDd::Zero();

    // precice::MeshHandle const& mesh =
    // _preciceInterface->getMeshHandle("BodyMesh");
    //
    VectorDs newPosition = _environmentForce * _dt * _dt
                           + 2.0 * _currentPosition -  _lastPosition;

    for (std::size_t i = 0; i < _vertexIds.size(); ++i) {
      for (unsigned d = 0; d < _dimensions; ++d) {
        _displacements[i * _dimensions + d]
          = newPosition(d) - _currentPosition(d);
      }
    }

    _lastPosition    = _currentPosition;
    _currentPosition = newPosition;

    _preciceInterface->writeBlockVectorData(_displacementsID,
                                            _vertexIds.size(),
                                            _vertexIds.data(),
                                            _displacements.data());
    logInfo("Write data {1}", _vertexIds.size());

    _dt = _preciceInterface->advance(_dt);

    return true;
  }

  Uni_PublicProperty(int,      type)
  Uni_PublicProperty(VectorDs, velocity)

private:
  unsigned                         _dimensions;
  Eigen::Matrix<long double, 3, 1> _environmentForce;
  Eigen::Matrix<long double, 3, 1> _lastPosition;
  Eigen::Matrix<long double, 3, 1> _currentPosition;
  double                           _dt;
  precice::SolverInterface*        _preciceInterface;
  std::vector<int>                 _vertexIds;
  std::vector<double>              _displacements;
  // VectorDs                             _shift;
  int _meshId;
  // int                                  _forcesID;
  int _displacementsID;
};
}
