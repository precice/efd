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

#include <map>
#include <unordered_map>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TGridGeometry,
          typename TMemory,
          typename TScalar,
          int TD,
          int TSolverType = 0>
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
      GridType,
      PpeStencilGenerator,
      PpeRhsGenerator,
      PpeResultAcquirer>
    PpeSolverType;

  typedef FluidSimulation::LinearSolver<
      GridType,
      VpeStencilGenerator,
      VpeRhsGenerator<0>,
      VpeResultAcquirer<0>>
    VxpeSolverType;

  typedef FluidSimulation::LinearSolver<
      GridType,
      VpeStencilGenerator,
      VpeRhsGenerator<1>,
      VpeResultAcquirer<1>>
    VypeSolverType;

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
          return CellAccessorType(i, &_memory, &_grid.innerGrid,
                                  &_gridGeometry);
        });

    _grid.initialize(localSize, cellAccessorFactory);

    // logGridInitializationInfo(_grid);
    logParallelTopologyInfo(_parallelDistribution);
  }

  int
  isOutside(int const& position) const {
    namespace pc = precice::constants;

    if (position == pc::positionOutsideOfGeometry()) {
      return 1;
    } else if (position == pc::positionOnGeometry()) {
      return 0;
    } else if (position == pc::positionInsideOfGeometry()) {
      return 0;
    }

    return -1;
  }

  int
  isTheSamePosition(int const& position0,
                    int const& position1) const {
    auto const& result0 = isOutside(position0);
    auto const& result1 = isOutside(position1);

    if ((result0 == -1)
        || (result1 == -1)) {
      return -1;
    }

    if (result0 == result1) {
      return 1;
    } else {
      return 0;
    }
  }

  int
  computeCellDistance(CellAccessorType const& accessor,
                      int const&              distance) const {
    int position = accessor.currentCell()->position();

    bool result = true;

    for (int currentDistance = 1;
         currentDistance <= distance;
         ++currentDistance) {
      for (int d = 0; d < TD; ++d) {
        for (int d2 = 0; d2 < 2; ++d2) {
          auto index = accessor.currentIndex();

          if (d2 == 0) {
            index(d) -= currentDistance;
          } else {
            index(d) += currentDistance;
          }

          if ((index(d) < _grid.innerGrid.leftIndent() (d))
              || (index(d) >= _grid.innerGrid.innerLimit() (d))) {
            continue;
          }

          auto const position1 = accessor.absoluteCell(index)->position();

          if (isTheSamePosition(position, position1) == 0) {
            result = false;

            return currentDistance;
          } else {
            result = true;
          }
        }
      }
    }

    return -1;
  }

  void
  initialize(precice::SolverInterface* preciceInteface,
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix) {
    _preciceInterface = preciceInteface;
    _ppeSolver.initialize(&_grid,
                          &_parallelDistribution,
                          &_parameters,
                          &_ghostHandler.ppeStencilGeneratorStack,
                          &_ghostHandler.ppeRhsGeneratorStack,
                          &_ghostHandler.ppeRhsAcquiererStack,
                          &_dt);

    _vxpeSolver.initialize(&_grid,
                           &_parallelDistribution,
                           &_parameters,
                           &_ghostHandler.vxpeStencilGeneratorStack,
                           &_ghostHandler.vxpeRhsGeneratorStack,
                           &_ghostHandler.vxpeRhsAcquiererStack,
                           &_dt);

    _vypeSolver.initialize(&_grid,
                           &_parallelDistribution,
                           &_parameters,
                           &_ghostHandler.vypeStencilGeneratorStack,
                           &_ghostHandler.vypeRhsGeneratorStack,
                           &_ghostHandler.vypeRhsAcquiererStack,
                           &_dt);

    _preciceInterface->initialize();
    _preciceInterface->initializeData();

    auto fluidMeshId = _preciceInterface->getMeshID("FluidMesh");

    for (auto const& accessor : _grid.innerGrid) {
      ImmersedBoundary::Fadlun::template
      computeDistances<CellAccessorType, TD>(accessor, _preciceInterface);
    }

    for (auto const& accessor : _grid.innerGrid) {
      Eigen::Vector4i pattern({ 1, 1, 1, 1 });

      bool doAdd = false;

      int distance = computeCellDistance(accessor, 2);
      // logInfo("{1} {2}", accessor.indexValues().transpose(),
      // distance);

      if (distance != -1) {
        for (int i = 0; i < 4; ++i) {
          auto const& position = pattern(i);

          if (!((i == distance)
                && (position == 1))) {
            continue;
          }

          if (i < 2) {
            if (isOutside(accessor.currentCell()->position()) == 1) {
              doAdd = true;
            }
          } else {
            if (isOutside(accessor.currentCell()->position()) == 0) {
              doAdd = true;
            }
          }
          break;
        }
      }

      if (doAdd) {
        VectorDsType position = accessor.currentPosition();
        position += 0.5 * accessor.currentWidth();
        auto vertexId =
          _preciceInterface->setMeshVertex(
            fluidMeshId,
            position.data());

        _vertexIds.insert(std::make_pair(accessor.indexValues(), vertexId));
      }
    }

    for (auto const& accessor : _grid) {
      accessor.currentCell()->velocity() = VelocityType::Zero();
      accessor.currentCell()->fgh()      = VelocityType::Zero();
      accessor.currentCell()->pressure() = 0.0;
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
                                                   _gridGeometry.
                                                   minCellWidth(),
                                                   _maxVelocity);
    _maxVelocity = VelocityType::Zero();

    if (_parallelDistribution.rank == 0) {
      logInfo("N = {1}; t = {2}", _iterationCount, _time);
    }

    auto fluidMeshId     = _preciceInterface->getMeshID("FluidMesh");
    auto fluidVertexSize = _preciceInterface->getMeshVertexSize(
      fluidMeshId);
    auto bodyMeshId     = _preciceInterface->getMeshID("Body");
    auto bodyVertexSize = _preciceInterface->getMeshVertexSize(bodyMeshId);

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

    int index = 0;

    for (auto const& accessor : _grid.innerGrid) {
      auto convection = ConvectionProcessing<TD>::compute(
        accessor,
        _parameters);

      auto previousConvection = accessor.currentCell()->convection();

      if (_iterationCount == 0) {
        previousConvection = convection;
      }

      auto diffusion = DiffusionProcessing<TD>::compute(accessor);

      diffusion = (1.0 / _parameters.re()) * diffusion;

      auto pressureGradient = computePressureGradient(accessor);

      VelocityType temp = VelocityType::Zero();

      VelocityType velocity;

      for (int d = 0; d < TD; ++d) {
        velocity(d) = 0.5 * (accessor.leftCellInDimension(d)->velocity(d)
                             + accessor.currentCell()->velocity(d));
      }

      if (_vertexIds.find(accessor.currentIndex()) != _vertexIds.end()) {
        temp = velocity + _dt * (
          -1.5 * convection + 0.5 * previousConvection
          + diffusion
          - pressureGradient);
        _preciceInterface->writeVectorData(
          fluidMeshVelocitiesId,
          _vertexIds[accessor.currentIndex()],
          temp.data());
      }

      accessor.currentCell()->convection() = convection;

      if (TSolverType == 1) {
        accessor.currentCell()->fgh()
          = accessor.currentCell()->velocity()
            + _dt * (-1.5 * convection + 0.5 * previousConvection
                     + 0.5 * diffusion
                     - pressureGradient
                     + _parameters.g());
      } else {
        accessor.currentCell()->fgh()
          = accessor.currentCell()->velocity()
            + _dt * (diffusion - convection + _parameters.g());
      }
    }

    _preciceInterface->advance(_dt);

    for (auto const& accessor : _grid.innerGrid) {
      auto find_it =  _vertexIds.find(accessor.currentIndex());

      if (find_it == _vertexIds.end()) {
        continue;
      }
      VelocityType force;
      _preciceInterface->readVectorData(
        fluidMeshForcesId,
        find_it->second,
        force.data());

      accessor.currentCell()->fgh() += force;
    }

    if (TSolverType == 1) {
      _vxpeSolver.solve();
      _vypeSolver.solve();
    } else {
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

  struct Comparison {
    bool
    operator()(VectorDiType const& one,
               VectorDiType const& two) const {
      return (one.array() == two.array()).all();
    }
  };

  struct Hash {
    std::size_t
    operator()(VectorDiType const& one) const {
      return one.sum();
    }
  };

  MemoryType                                  _memory;
  GridGeometryType                            _gridGeometry;
  GridType                                    _grid;
  std::unordered_map<VectorDiType, int, Hash> _vertexIds;
  ParametersType                              _parameters;
  ParallelDistributionType                    _parallelDistribution;
  precice::SolverInterface*                   _preciceInterface;
  TScalar                                     _dt;
  TScalar                                     _time;
  TScalar                                     _timeLimit;
  TScalar                                     _lastPlotTimeStamp;
  TScalar                                     _plotInterval;
  int                                         _iterationCount;
  int                                         _iterationLimit;
  PpeSolverType                               _ppeSolver;
  VxpeSolverType                              _vxpeSolver;
  VypeSolverType                              _vypeSolver;
  VelocityType                                _maxVelocity;
  GhostHandlersType                           _ghostHandler;
  VtkPlotType                                 _plot;
};
}
}
#endif
