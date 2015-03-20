#pragma once

#include "Configuration.hpp"
#include "GhostLayer/SfsfdHandlers.hpp"
#include "ImmersedBoundary/BodyForce/functions.hpp"
#include "ImmersedBoundary/geometryfunctions.hpp"
#include "LinearSolver.hpp"
#include "Private/mpigenerics.hpp"
#include "Private/unorderedmapgenerics.hpp"
#include "Solver.hpp"
#include "SolverTraits.hpp"
#include "functions.hpp"

#include <Uni/Logging/macros>

#include <map>
#include <unordered_map>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TGridGeometry,
          typename TSimulationResultsWriter,
          typename TScalar,
          int TD>
class SfsfdSolver : public Solver {
public:
  using BaseType = Solver;

  using Path = typename BaseType::Path;

  using SolverTraitsType = SolverTraits<TGridGeometry,
                                        TScalar,
                                        TD>;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridGeometryType = typename SolverTraitsType::GridGeometryType;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

  using GhostHandlersType
          = typename GhostLayer::SfsfdHandlers<Dimensions>;

  typedef PpeStencilGenerator<CellAccessorType, ParallelDistributionType>
    PpeStencilGeneratorType;
  typedef PpeRhsGenerator<CellAccessorType>
    PpeRhsGeneratorType;
  typedef PpeResultAcquirer1<CellAccessorType>
    PpeResultAcquirerType;

  typedef FluidSimulation::LinearSolver<
      GridType,
      PpeStencilGeneratorType,
      PpeRhsGeneratorType,
      PpeResultAcquirerType>
    PpeSolverType;

public:
  SfsfdSolver() {}

  SfsfdSolver(SfsfdSolver const& other) = delete;

  ~SfsfdSolver() {}

  SfsfdSolver const&
  operator=(SfsfdSolver const& other) = delete;

  void
  initialize(precice::SolverInterface* preciceInteface,
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix) {
    _preciceInterface = preciceInteface;

    _ppeStencilGenerator.initialize(_memory.parallelDistribution());
    _ppeRhsGenerator.initialize(&_memory.timeStepSize());

    _ppeSolver.initialize(_memory.grid(),
                          _memory.parallelDistribution(),
                          &_ppeStencilGenerator,
                          &_ppeRhsGenerator,
                          &_ppeResultAcquierer,
                          &_ghostHandlers.ppeStencilGeneratorStack,
                          &_ghostHandlers.ppeRhsGeneratorStack,
                          &_ghostHandlers.ppeRhsAcquiererStack);

    _resultWriter.initialize(&_memory, outputDirectory, fileNamePrefix);

    _preciceInterface->initialize();
    _preciceInterface->initializeData();

    ImmersedBoundary::compute_position_in_respect_to_geometry(_memory.grid(),
                                                              _preciceInterface);

    auto fluidMeshId = _preciceInterface->getMeshID("FluidMesh");

    for (auto const& accessor : _memory.grid()->innerGrid) {
      // logInfo("{1} {2}", accessor.index().transpose(),
      //         accessor.globalIndex());

      int distance
        = ImmersedBoundary::compute_cell_layer_along_geometry_interface(
        accessor, 2);

      bool doAdd
        = ImmersedBoundary::validate_layer_number(
        accessor, distance, 2, 2);

      if (doAdd) {
        VectorDsType position = accessor.pressurePosition();
        auto         vertexId
          = _preciceInterface->setMeshVertex(fluidMeshId, position.data());

        _vertexIds.insert(std::make_pair(accessor.indexValues(), vertexId));
      }
    }

    for (auto& accessor : * _memory.grid()) {
      accessor.velocity() = VectorDsType::Zero();
      accessor.fgh()      = VectorDsType::Zero();
      accessor.pressure() = 0.0;
    }

    _ghostHandlers.executeVelocityMpiExchange();

    _lastPlotTimeStamp = 0.0;

    _memory.maxVelocity()     = VectorDsType::Zero();
    _memory.timeStepSize()    = 1.0;
    _memory.time()            = 0.0;
    _memory.iterationNumber() = 0;

    if (_plotInterval >= 0) {
      _resultWriter.writeGeometry();
      _resultWriter.writeAttributes();
    }
  }

  bool
  iterate() {
    if (_timeLimit > 0 && _memory.time() >= _timeLimit) {
      return false;
    }

    if (_iterationLimit > 0 && _memory.iterationNumber() >= _iterationLimit) {
      return false;
    }

    _memory.timeStepSize() = TimeStepProcessing<TScalar, TD>::compute(
      _memory.parameters()->re(),
      _memory.parameters()->tau(),
      _memory.gridGeometry()->minCellWidth(),
      _memory.maxVelocity());
    _memory.maxVelocity() = VectorDsType::Zero();

    if (_memory.parallelDistribution()->rank == 0) {
      logInfo("N = {1}; t = {2}", _memory.iterationNumber(), _memory.time());
    }

    auto fluidMeshId
      = _preciceInterface->getMeshID("FluidMesh");
    auto fluidVertexSize
      = _preciceInterface->getMeshVertexSize(fluidMeshId);
    auto bodyMeshId
      = _preciceInterface->getMeshID("Body");
    auto bodyVertexSize
      = _preciceInterface->getMeshVertexSize(bodyMeshId);
    auto fluidMeshVelocitiesId
      = _preciceInterface->getDataID("Velocities",
                                     _preciceInterface->getMeshID("FluidMesh"));
    auto fluidMeshForcesId
      = _preciceInterface->getDataID("Forces",
                                     _preciceInterface->getMeshID("FluidMesh"));
    auto bodyMeshVelocitiesId
      = _preciceInterface->getDataID("Velocities",
                                     _preciceInterface->getMeshID("Body"));
    auto bodyMeshForcesId
      = _preciceInterface->getDataID("Forces",
                                     _preciceInterface->getMeshID("Body"));

    int index = 0;

    for (auto const& accessor : _memory.grid()->innerGrid) {
      auto convection = ConvectionProcessing<TD>::compute(
        accessor,
        _memory.parameters());

      auto previousConvection = accessor.convection();

      if (_memory.iterationNumber() == 0) {
        previousConvection = convection;
      }

      accessor.convection() = convection;

      auto diffusion = DiffusionProcessing<TD>::compute(accessor);

      diffusion = (1.0 / _memory.parameters()->re()) * diffusion;

      VectorDsType temp = VectorDsType::Zero();

      VectorDsType velocity;

      for (int d = 0; d < TD; ++d) {
        velocity(d) = 0.5 * (accessor.velocity(d, -1, d) +
                             accessor.velocity(d));
      }

      auto pressureGradient = computePressureGradient(accessor);

      if (_vertexIds.find(accessor.index()) != _vertexIds.end()) {
        temp = velocity + _memory.timeStepSize() * (
          -1.5 * convection + 0.5 * previousConvection
          + diffusion
          - pressureGradient);

        _preciceInterface->writeVectorData(
          fluidMeshVelocitiesId,
          _vertexIds[accessor.index()],
          temp.data());
      }

      accessor.fgh()
        = accessor.velocity()
          + _memory.timeStepSize() * (diffusion - convection +
                                      _memory.parameters()->g());
    }

    _preciceInterface->advance(_memory.timeStepSize());

    for (auto& accessor : _memory.grid()->innerGrid) {
      auto find_it =  _vertexIds.find(accessor.index());

      if (find_it == _vertexIds.end()) {
        continue;
      }
      VectorDsType force;
      _preciceInterface->readVectorData(fluidMeshForcesId,
                                        find_it->second,
                                        force.data());

      accessor.fgh() += force;
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandlers.fghInitialization[d][d2]();
      }
    }

    _ghostHandlers.executeFghMpiExchange();

    _ppeSolver.solve();

    _ghostHandlers.executePressureMpiExchange();

    for (auto& accessor : _memory.grid()->innerGrid) {
      for (int d = 0; d < TD; ++d) {
        accessor.velocity(d)
          = accessor.fgh(d)
            - _memory.timeStepSize() /
            (0.5 * (accessor.width(d, +1, d)
                    + accessor.width() (d)))
            * (accessor.pressure(d, +1)
               - accessor.pressure());
      }
      computeMaxVelocity<CellAccessorType, TScalar, TD>
        (accessor, _memory.maxVelocity());
    }

    _ghostHandlers.executeVelocityInitialization();

    _ghostHandlers.executeVelocityMpiExchange();

    VectorDsType force = VectorDsType::Zero();

    for (auto& accessor : _memory.grid()->innerGrid) {
      ImmersedBoundary::BodyForce::
      template computeCellForce<CellAccessorType>(
        accessor,
        _memory.parameters()->re(),
        force);
    }

    Private::mpiAllReduce<TScalar>(MPI_IN_PLACE,
                                   force.data(),
                                   TD,
                                   MPI_SUM,
                                   PETSC_COMM_WORLD);

    _memory.time() += _memory.timeStepSize();
    ++_memory.iterationNumber();

    if (_plotInterval >= 0) {
      if ((_memory.time() - _lastPlotTimeStamp) > _plotInterval) {
        _lastPlotTimeStamp = _memory.time();

        _resultWriter.writeAttributes();
      }
    }

    force *= 2 / (0.3 * 0.3 * 0.1);

    if (_memory.parallelDistribution()->rank == 0) {
      logInfo("{1}", force.transpose());
    }

    return true;
  }

  struct Comparison {
    bool
    operator()(VectorDiType const& one,
               VectorDiType const& two) const {
      return (one.array() == two.array()).all();
    }
  };

  MemoryType _memory;
  std::unordered_map < VectorDiType, int, Private::Hash < VectorDiType >>
  _vertexIds;
  precice::SolverInterface* _preciceInterface;
  long double               _timeLimit;
  TScalar                   _lastPlotTimeStamp;
  TScalar                   _plotInterval;
  unsigned long long        _iterationLimit;

  PpeStencilGeneratorType  _ppeStencilGenerator;
  PpeRhsGeneratorType      _ppeRhsGenerator;
  PpeResultAcquirerType    _ppeResultAcquierer;
  PpeSolverType            _ppeSolver;
  GhostHandlersType        _ghostHandlers;
  TSimulationResultsWriter _resultWriter;
};
}
}
