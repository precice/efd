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
    _type          = configuration->get<unsigned>("/Type");
    _isPreciceMode = configuration->get<bool>("/PreciceMode");
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

    unsigned vertices_size = _preciceInterface->getMeshVertexSize(_meshId);

    _vertexIds.resize(vertices_size);
    _displacements.resize(_dimensions * vertices_size);
    _displacementsID = _preciceInterface->getDataID("Displacements", _meshId);

    if (_isPreciceMode) {
      logInfo("Precice mode is on");
      _forcesID = _preciceInterface->getDataID("Forces", _meshId);
      _forces.resize(_dimensions * vertices_size);
    }

    for (unsigned i = 0; i < vertices_size; ++i) {
      _vertexIds[i] = i;
    }
    _lastPosition    = VectorDs::Zero();
    _currentPosition = VectorDs::Zero();

    if (_type == 1) {
      _velocity       = _environmentForce * 0.1;
      _type1Direction = true;
    }
  }

  bool
  iterate() {
    if (!_preciceInterface->isCouplingOngoing()) {
      return false;
    }

    if (_type == 0) {
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
      // logInfo("Write data {1}", _vertexIds.size());
    } else if (_type == 1) {
      // const double PI = 3.141592653589793238463;
      VectorDs     newPosition;

      for (unsigned d = 0; d < _dimensions; ++d) {
        _velocity(d)
          = _environmentForce(d);// * 0.9
            // * std::sin(PI * _currentPosition(d) / 1.7) + _environmentForce(d) * 0.1;

        if (_type1Direction) {
          newPosition(d) = _currentPosition(d) + _velocity(d) * _dt;
        } else {
          newPosition(d) = _currentPosition(d) - _velocity(d) * _dt;
        }

        if (newPosition(d) >= 1.7) {
          _type1Direction = false;
        }

        if (newPosition(d) <= 0.0) {
          _type1Direction = true;
        }

        if (_type1Direction) {
          newPosition(d) = _currentPosition(d) + _velocity(d) * _dt;
        } else {
          newPosition(d) = _currentPosition(d) - _velocity(d) * _dt;
        }
      }

      for (std::size_t i = 0; i < _vertexIds.size(); ++i) {
        for (unsigned d = 0; d < _dimensions; ++d) {
          _displacements[i * _dimensions + d] = newPosition(d) - _currentPosition(d);

          if (_isPreciceMode) {
            _forces[i * _dimensions + d] = _velocity(d);
          }
        }
      }
      _currentPosition = newPosition;

      _preciceInterface->writeBlockVectorData(_displacementsID,
                                              _vertexIds.size(),
                                              _vertexIds.data(),
                                              _displacements.data());

      if (_isPreciceMode) {
        logInfo("Write forces");
        _preciceInterface->writeBlockVectorData(_forcesID,
                                                _vertexIds.size(),
                                                _vertexIds.data(),
                                                _forces.data());
      }
    }

    _dt = _preciceInterface->advance(_dt);

    return true;
  }

  Uni_PublicProperty(unsigned, type)
  Uni_PublicProperty(bool,     isPreciceMode)
  Uni_PublicProperty(VectorDs, velocity)

private:
  unsigned                         _dimensions;
  Eigen::Matrix<long double, 3, 1> _environmentForce;
  Eigen::Matrix<long double, 3, 1> _lastPosition;
  Eigen::Matrix<long double, 3, 1> _currentPosition;
  double                           _dt;
  bool                             _type1Direction;
  precice::SolverInterface*        _preciceInterface;
  std::vector<int>                 _vertexIds;
  std::vector<double>              _displacements;
  std::vector<double>              _forces;
  int                              _meshId;
  int                              _displacementsID;
  int                              _forcesID;
};
}
