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
  using Scalar = long double;

  using VectorDd = Eigen::Matrix<double, 3, 1>;

  using VectorDs = Eigen::Matrix<Scalar, 3, 1>;

public:
  Simulation(Configuration const* configuration) {
    _dimensions = configuration->get<unsigned>("/Dimensions");

    if (!configuration->isOfType<bool>("/EnvironmentForce")) {
      _environmentForce = configuration->get<VectorDs>("/EnvironmentForce");
    } else {
      _environmentForce = VectorDs::Zero();
    }
    _type          = configuration->get<unsigned>("/Type");
    _mass          = configuration->get<Scalar>("/Mass");
    _isPreciceMode = configuration->get<bool>("/PreciceMode");
    _positionLimit = configuration->get<VectorDs>("/PositionLimit");
  }

  Simulation(Simulation const&) = delete;

  ~Simulation() {}

  Simulation&
  operator=(Simulation const&) = delete;

  void
  initialize(precice::SolverInterface* precice_interface) {
    // _type = type;
    _preciceInterface = precice_interface;
    // logInfo("PreCICE's initialize methods is being invoked ...");
    _dt = _preciceInterface->initialize();
    // logInfo("PreCICE's initialize methods has been finished");
    _dt = 0.0;

    if (_dimensions == 0u) {
      _dimensions = _preciceInterface->getDimensions();
    }

    if (static_cast<int>(_dimensions) != _preciceInterface->getDimensions()) {
      throwException("Inconsistent dimension size of the problem: "
                     "structure's configuration dimensions ('{1}') is not equal to "
                     "PreCICE's configuration dimensions ('{2}')",
                     _dimensions,
                     _preciceInterface->getDimensions());
    }

    if (!_preciceInterface->hasMesh("BodyMesh")) {
      throwException("Precice does not have 'BodyMesh' in its configuration");
    }

    _ibMeshId = _preciceInterface->getMeshID("BodyMesh");

    if (_isPreciceMode) {
      if (!_preciceInterface->hasMesh("CouplingBodyMesh")) {
        throwException("Precice does not have 'CouplingBodyMesh' in its configuration");
      }
      _couplingMeshId = _preciceInterface->getMeshID("CouplingBodyMesh");
    } else {
      _couplingMeshId = _preciceInterface->getMeshID("BodyMesh");
    }

    unsigned vertices_size = _preciceInterface->getMeshVertexSize(_ibMeshId);

    _vertexIds.resize(vertices_size);
    _displacements.resize(_dimensions * vertices_size);
    _displacementsID = _preciceInterface->getDataID("Displacements", _ibMeshId);

    if (_isPreciceMode) {
      // logInfo("Precice mode is on");
      // _forcesID = _preciceInterface->getDataID("Forces", _ibMeshId);
      // _forces.resize(_dimensions * vertices_size);
    }

    for (unsigned i = 0; i < vertices_size; ++i) {
      _vertexIds[i] = i;
    }
    _previousPosition = VectorDs::Zero();
    _currentPosition  = VectorDs::Zero();

    if (_type == 1) {
      _velocity       = _environmentForce * 0.1;
      _type1Direction = true;
    }
  }

  bool
  iterate() {
    namespace pc = precice::constants;
    static std::string const writeCheckpoint(pc::actionWriteIterationCheckpoint());
    static std::string const readCheckpoint(pc::actionReadIterationCheckpoint());

    if (!_preciceInterface->isCouplingOngoing()) {
      return false;
    }

    while (_preciceInterface->isCouplingOngoing()) {
      VectorDs previous_position;
      VectorDs current_position;
      bool     type_1_direction;

      if (_preciceInterface->isActionRequired(writeCheckpoint)) {
        previous_position = _previousPosition;
        current_position  = _currentPosition;
        type_1_direction  = _type1Direction;
        _preciceInterface->fulfilledAction(writeCheckpoint);
      }

      // logInfo("Start of iteration");

      if (_type == 0) {
        if ((_currentPosition.cwiseAbs().array()
             > _positionLimit.cwiseAbs().array()).any()) {
          return false;
        }

        if (!_preciceInterface->hasData("CouplingForces", _couplingMeshId)) {
          throwException("Precice configuration does not have 'CouplingForces' data"
                         " related to 'BodyMesh'");
        }
        auto const stressesId = _preciceInterface->getDataID("CouplingForces",
                                                             _couplingMeshId);

        unsigned vertices_size = _preciceInterface->getMeshVertexSize(_couplingMeshId);

        std::vector<double> stresses;
        stresses.resize(_dimensions * vertices_size);

        _preciceInterface->readBlockVectorData(stressesId,
                                               _vertexIds.size(),
                                               _vertexIds.data(),
                                               stresses.data());

        VectorDs force = _environmentForce;

        for (std::size_t i = 0; i < _vertexIds.size(); ++i) {
          for (unsigned d = 0; d < _dimensions; ++d) {
            // logInfo("Stresses = {1}", stresses[i * _dimensions + d]);
            force(d) += stresses[i * _dimensions + d];
          }
        }

        VectorDs newPosition = force / _mass * _dt * _dt
                               + 2.0 * _currentPosition - _previousPosition;

        logInfo("Force = {1} {2}", force.transpose(), _mass);
        logInfo("Position = {1}",  newPosition.transpose());

        for (std::size_t i = 0; i < _vertexIds.size(); ++i) {
          for (unsigned d = 0; d < _dimensions; ++d) {
            _displacements[i * _dimensions + d] = newPosition(d) - _currentPosition(d);
          }
        }

        _previousPosition = _currentPosition;
        _currentPosition  = newPosition;

        _preciceInterface->writeBlockVectorData(_displacementsID,
                                                _vertexIds.size(),
                                                _vertexIds.data(),
                                                _displacements.data());

        if (_isPreciceMode) {
          if (!_preciceInterface->hasData("Displacements", _couplingMeshId)) {
            throwException("Precice configuration does not have 'Displacements' data"
                           " related to 'CouplingBodyMesh'");
          }
          auto const coupling_mesh_displacements_id
            = _preciceInterface->getDataID("Displacements", _couplingMeshId);
          _preciceInterface->writeBlockVectorData(coupling_mesh_displacements_id,
                                                  _vertexIds.size(),
                                                  _vertexIds.data(),
                                                  _displacements.data());
        }
      } else if (_type == 1) {
        // const double PI = 3.141592653589793238463;
        VectorDs newPosition;

        for (unsigned d = 0; d < _dimensions; ++d) {
          _velocity(d)
            = _environmentForce(d); // * 0.9
          // * std::sin(PI * _currentPosition(d) / 1.7) + _environmentForce(d) *
          // 0.1;

          if (!_type1Direction) {
            _velocity(d) *= -1;
          }

          newPosition(d) = _currentPosition(d) + _velocity(d) * _dt;

          if (newPosition(d) > _positionLimit(d)) {
            _type1Direction = false;
            _velocity(d)   *= -1;
            newPosition(d)  = _currentPosition(d) + _velocity(d) * _dt;
          }

          if (newPosition(d) < 0.0) {
            _type1Direction = true;
            _velocity(d)   *= -1;
            newPosition(d)  = _currentPosition(d) + _velocity(d) * _dt;
          }
        }

        // logInfo("{1} {2} {3}", newPosition(0), _currentPosition(0),
        // _velocity(0));

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

        // if (_isPreciceMode) {
        // logInfo("Write forces {1}", _velocity.transpose());
        // _preciceInterface->writeBlockVectorData(_forcesID,
        // _vertexIds.size(),
        // _vertexIds.data(),
        // _forces.data());
        // }
      } else if (_type == 2) {
        VectorDs newPosition = _environmentForce * _dt * _dt
                               + 2.0 * _currentPosition -  _previousPosition;

        for (std::size_t i = 0; i < _vertexIds.size(); ++i) {
          for (unsigned d = 0; d < _dimensions; ++d) {
            _displacements[i * _dimensions + d]
              = newPosition(d) - _currentPosition(d);
          }
        }

        _previousPosition = _currentPosition;
        _currentPosition  = newPosition;

        _preciceInterface->writeBlockVectorData(_displacementsID,
                                                _vertexIds.size(),
                                                _vertexIds.data(),
                                                _displacements.data());
        // logInfo("Write data {1}", _vertexIds.size());
      }

      // logInfo("PreCICE's advance methods is being invoked ...");
      _dt = _preciceInterface->advance(_dt);
      // logInfo("PreCICE's advance methods has been finished");

      if (_preciceInterface->isActionRequired(readCheckpoint)) {
        _previousPosition = previous_position;
        _currentPosition  = current_position;
        _type1Direction   = type_1_direction;
        _preciceInterface->fulfilledAction(readCheckpoint);
      } else {
        break;
      }

      // logInfo("End of iteration");
    }

    return true;
  }

  Uni_PublicProperty(unsigned, type)
  Uni_PublicProperty(bool,     isPreciceMode)
  Uni_PublicProperty(VectorDs, velocity)

private:
  unsigned                  _dimensions;
  VectorDs                  _environmentForce;
  VectorDs                  _positionLimit;
  VectorDs                  _previousPosition;
  VectorDs                  _currentPosition;
  Scalar                    _mass;
  Scalar                    _dt;
  bool                      _type1Direction;
  precice::SolverInterface* _preciceInterface;
  std::vector<int>          _vertexIds;
  std::vector<double>       _displacements;
  std::vector<double>       _forces;
  int                       _ibMeshId;
  int                       _couplingMeshId;
  int                       _displacementsID;
  // int                       _forcesID;
};
}
