#ifndef FsiSimulation_FluidSimulation_FdSimulation_hpp
#define FsiSimulation_FluidSimulation_FdSimulation_hpp

#include "Configuration.hpp"
#include "GhostLayer/Handlers.hpp"
#include "Grid.hpp"
#include "ImmersedBoundary/BodyForce/functions.hpp"
#include "ImmersedBoundary/Fadlun/functions.hpp"
#include "ImmersedBoundary/FeedbackForcing/functions.hpp"
#include "LinearSolver.hpp"
#include "ParallelDistribution.hpp"
#include "Parameters.hpp"
#include "Private/mpigenerics.hpp"
#include "Simulation.hpp"
#include "VtkPlot.hpp"
#include "functions.hpp"

#include <Uni/Logging/macros>
#include <Uni/Stopwatch>

#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TGridGeometry,
          typename TMemory,
          typename TScalar,
          int TD>
class FdSimulation : public Simulation {
public:
  typedef Simulation              BaseType;
  typedef typename BaseType::Path Path;

  typedef TScalar                                     ScalarType;
  typedef Grid<TMemory, TGridGeometry, TD>            GridType;
  typedef typename GridType::CellAccessorType         CellAccessorType;
  typedef typename CellAccessorType::MemoryType       MemoryType;
  typedef typename CellAccessorType::CellType         CellType;
  typedef typename CellAccessorType::GridGeometryType GridGeometryType;
  typedef typename CellType::Velocity                 VelocityType;
  typedef typename GridType::VectorDi                 VectorDiType;
  typedef typename GridGeometryType::VectorDs         VectorDsType;
  typedef typename GridType::Factory                  CellAccessorFactory;

  typedef Parameters<TScalar, TD> ParametersType;

  typedef ParallelDistribution<TD> ParallelDistributionType;

  typedef FluidSimulation::LinearSolver<
      TMemory,
      TGridGeometry,
      PpeStencilGenerator<TD>,
      ScalarType,
      TD>
    PpeSolverType;

  typedef typename GhostLayer::Handlers<TD> GhostHandlersType;

  typedef
    VtkPlot<TMemory, TGridGeometry, TScalar, TD>
    VtkPlotType;

public:
  FdSimulation() {}

  FdSimulation(FdSimulation const& other) = delete;

  ~FdSimulation() {}

  FdSimulation const&
  operator=(FdSimulation const& other) = delete;

  void
  setParameters(VectorDiType const& localSize,
                VectorDsType const& width) {
    _gridGeometry.initialize(
      width,
      _parallelDistribution.globalCellSize,
      _parallelDistribution.corner);

    _memory.allocate(localSize);

    CellAccessorFactory cellAccessorFactory(
      [&] (VectorDiType const& i) {
        return CellAccessorType(i, &_memory, &_grid.innerGrid, &_gridGeometry);
      });

    _grid.initialize(localSize, cellAccessorFactory);

    // logGridInitializationInfo(_grid);
    logParallelTopologyInfo(_parallelDistribution);
  }

  void
  initialize(precice::SolverInterface* preciceInteface,
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix) {
    _preciceInterface = preciceInteface;
    _ppeSolver.initialize(&_grid,
                          &_parallelDistribution,
                          &_ghostHandler,
                          &_dt);

    auto fluidMeshId = _preciceInterface->getMeshID("FluidMesh");

    for (auto const& accessor : _grid) {
      accessor.currentCell()->velocity() = VelocityType::Zero();
      accessor.currentCell()->fgh()      = VelocityType::Zero();
      accessor.currentCell()->pressure() = 0.0;

      bool isInnerCell = true;

      for (int d = 0; d < TD; ++d) {
        if ((accessor.indexValue(d) >= _grid.innerGrid.innerLimit(d))
            || accessor.indexValue(d) < _grid.innerGrid.leftIndent(d)) {
          isInnerCell = false;
          break;
        }
      }

      if (isInnerCell) {
        VectorDsType position = accessor.currentPosition();
        position += 0.5 * accessor.currentWidth();
        auto vertexId =
          _preciceInterface->setMeshVertex(
            fluidMeshId,
            position.data());
        _vertexIds.push_back(vertexId);
      }

      // ImmersedBoundary::Fadlun::template computeDistances<CellAccessorType,
      // TD>(
      // accessor,
      // _preciceInterface);
    }
        logInfo("@@@@@@@@@@@@@");

    _preciceInterface->initialize();
    _preciceInterface->initializeData();

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.velocityInitialization[d][d2]();
      }
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.mpiVelocityExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _maxVelocity       = VelocityType::Zero();
    _dt                = 1.0;
    _time              = 0.0;
    _lastPlotTimeStamp = 0.0;
    _iterationCount    = 0;

    if (_plotInterval >= 0) {
      _plot.initialize(&_grid,
                       &_parallelDistribution,
                       &_gridGeometry,
                       outputDirectory,
                       fileNamePrefix);

      _plot.plot(_iterationCount, _time, _dt);
    }
  }

  bool
  iterate() {
    if (_timeLimit > 0 && _time >= _timeLimit) {
      return false;
    }

    if (_iterationLimit > 0 && _iterationCount >= _iterationLimit) {
      return false;
    }

    _dt = TimeStepProcessing<TScalar, TD>::compute(_parameters.re(),
                                                   _parameters.tau(),
                                                   _gridGeometry.minCellWidth(),
                                                   _maxVelocity);
    _maxVelocity = VelocityType::Zero();

    if (_parallelDistribution.rank == 0) {
        logInfo("N = {1}; t = {2}", _iterationCount, _time);
    }
    // logInfo("Time step size {1}",   _dt);

    // for (auto accessor : _grid.innerGrid) {
    // for (int d = 0; d < TD; ++d) {
    // if (!ImmersedBoundary::Fadlun::treatBoundary
    // <CellAccessorType, TScalar, TD>(
    // accessor,
    // _preciceInterface,
    // _dt,
    // d)) {
    ////
    // } else {
    // computeMaxVelocity<CellAccessorType, TScalar, TD>
    // (accessor, _maxVelocity);
    // }
    // }
    // }

    // for (auto const& accessor : _grid.innerGrid) {
    // _preciceInterface->mapWriteDataFrom(
    // _preciceInterface->getMeshID("MeshName"));
    // }
    //
    auto fluidMeshId     = _preciceInterface->getMeshID("FluidMesh");
    auto fluidVertexSize = _preciceInterface->getMeshVertexSize(fluidMeshId);
    auto bodyMeshId      = _preciceInterface->getMeshID("Body");
    auto bodyVertexSize  = _preciceInterface->getMeshVertexSize(bodyMeshId);

    auto fluidMeshVelocitiesId = _preciceInterface->getDataID(
      "Velocities",
      _preciceInterface->getMeshID("FluidMesh"));
    auto fluidMeshForcesId = _preciceInterface->getDataID(
      "Forces",
      _preciceInterface->getMeshID("FluidMesh"));
    auto bodyMeshVelocitiesId = _preciceInterface->getDataID(
      "Velocities",
      _preciceInterface->getMeshID("Body"));
    auto bodyMeshForcesId = _preciceInterface->getDataID(
      "Forces",
      _preciceInterface->getMeshID("Body"));

    // logInfo("FluidMesh1 {1}", fluidMeshId);
    // logInfo("FluidMesh {1}",  bodyMeshId);
    // logInfo("Body Size{1}",   bodyVertexSize);

    int index = 0;

    for (auto const& accessor : _grid.innerGrid) {
      // typedef FghProcessing<TD> Fgh;
      // Fgh::compute(accessor, _parameters, _dt);

      accessor.currentCell()->fgh() = VelocityType::Zero();

      // bool isBody = false;

      // for (int d = 0; d < TD; ++d) {
      // if (ImmersedBoundary::FeedbackForcing
      // ::template treatBoundary<CellAccessorType, TD>(accessor,
      // _parameters.alpha(),
      // d)) {
      // isBody = true;
      // }
      // }

      auto convection = ConvectionProcessing<TD>::compute(
        accessor,
        _parameters);

      auto previousConvection = accessor.currentCell()->convection();

      if (_iterationCount == 0) {
        previousConvection = convection;
      }

      auto diffusion = DiffusionProcessing<TD>::compute(accessor);

      diffusion = (1.0 / _parameters.re()) * diffusion;
      VelocityType temp = VelocityType::Zero();

      // if (isBody) {
      temp = -accessor.currentCell()->velocity() + _dt * (
        1.5 * convection - 0.5 * previousConvection
        - diffusion
        + computePressureGradient(accessor));
      // }

      _preciceInterface->writeVectorData(
        fluidMeshVelocitiesId,
        _vertexIds[index],
        temp.data());

      accessor.currentCell()->convection() = convection;

      accessor.currentCell()->fgh()
        += accessor.currentCell()->velocity()
           + _dt * (diffusion - convection + _parameters.g());
      ++index;
    }

    _preciceInterface->mapWriteDataFrom(fluidMeshId);

    index = 0;

    for (; index < _preciceInterface->getMeshVertexSize(bodyMeshId);
         ++index) {
      VelocityType velocity;
      _preciceInterface->readVectorData(
        bodyMeshVelocitiesId,
        index,
        velocity.data());

      velocity = velocity / _dt;

      //if (velocity != VelocityType::Zero()) {
      //  logInfo("{1}", velocity.transpose());
      //}

      _preciceInterface->writeVectorData(
        bodyMeshForcesId,
        index,
        velocity.data());
    }

    _preciceInterface->mapReadDataTo(fluidMeshId);

    index = 0;

    for (auto const& accessor : _grid.innerGrid) {
      VelocityType velocity;
      _preciceInterface->readVectorData(
        fluidMeshForcesId,
        _vertexIds[index],
        velocity.data());

      if (velocity != VelocityType::Zero()) {
        logInfo("{2} {1}", velocity.transpose(), accessor.indexValues().transpose());
      }

      accessor.currentCell()->fgh() += _dt * velocity;

      ++index;
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.fghInitialization[d][d2]();
      }
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.mpiFghExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    _ppeSolver.solve();

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.mpiPressureExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    for (auto const& accessor : _grid.innerGrid) {
      typedef VelocityProcessing<CellAccessorType,
                                 TScalar,
                                 TD> VelocityProcessingType;

      for (int d = 0; d < TD; ++d) {
        VelocityProcessingType::compute(accessor, d, _dt);
      }
      computeMaxVelocity<CellAccessorType, TScalar, TD>
        (accessor, _maxVelocity);
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.velocityInitialization[d][d2]();
      }
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.mpiVelocityExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    VectorDsType force = VectorDsType::Zero();

    for (auto const& accessor : _grid.innerGrid) {
      ImmersedBoundary::BodyForce::
      template computeCellForce<CellAccessorType, TD>(
        accessor,
        _parameters.re(),
        force);
    }

    Private::mpiAllReduce<TScalar>(MPI_IN_PLACE,
                                   force.data(),
                                   TD,
                                   MPI_SUM,
                                   PETSC_COMM_WORLD);

    _time += _dt;
    ++_iterationCount;

    if (_plotInterval >= 0) {
      if ((_time - _lastPlotTimeStamp) > _plotInterval) {
        _lastPlotTimeStamp = _time;
        _plot.plot(_iterationCount, _time, _dt);
      }
    }

    force *= 2 / (0.3 * 0.3 * 0.1);

    if (_parallelDistribution.rank == 0) {
      logInfo("{1}", force.transpose());
    }

    return true;
  }

  MemoryType                _memory;
  GridGeometryType          _gridGeometry;
  GridType                  _grid;
  std::vector<int>          _vertexIds;
  ParametersType            _parameters;
  ParallelDistributionType  _parallelDistribution;
  precice::SolverInterface* _preciceInterface;
  TScalar                   _dt;
  TScalar                   _time;
  TScalar                   _timeLimit;
  TScalar                   _lastPlotTimeStamp;
  TScalar                   _plotInterval;
  int                       _iterationCount;
  int                       _iterationLimit;
  PpeSolverType             _ppeSolver;
  VelocityType              _maxVelocity;
  GhostHandlersType         _ghostHandler;
  VtkPlotType               _plot;
};
}
}
#endif
