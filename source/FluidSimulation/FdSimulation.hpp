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
#include "XdmfHdf5Output/XdmfHdf5Writer.hpp"
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

  typedef PpeStencilGenerator<CellAccessorType, ParallelDistributionType>
    PpeStencilGeneratorType;
  typedef PpeRhsGenerator<CellAccessorType>
    PpeRhsGeneratorType;
  typedef PpeResultAcquirer1<CellAccessorType>
    PpeResultAcquirer1Type;
  typedef PpeResultAcquirer2<CellAccessorType>
    PpeResultAcquirer2Type;

  typedef VpeStencilGenerator
    <CellAccessorType, ParallelDistributionType, ParametersType>
    VpeStencilGeneratorType;

  typedef VpeRhsGenerator<CellAccessorType, 0>
    VxpeRhsGeneratorType;
  typedef VpeResultAcquirer<0>
    VxpeResultAcquirerType;

  typedef VpeRhsGenerator<CellAccessorType, 1>
    VypeRhsGeneratorType;
  typedef VpeResultAcquirer<1>
    VypeResultAcquirerType;

  typedef FluidSimulation::LinearSolver<
      GridType,
      PpeStencilGeneratorType,
      PpeRhsGeneratorType,
      PpeResultAcquirer1Type>
    PpeSolver1Type;

  typedef FluidSimulation::LinearSolver<
      GridType,
      PpeStencilGeneratorType,
      PpeRhsGeneratorType,
      PpeResultAcquirer2Type>
    PpeSolver2Type;

  typedef FluidSimulation::LinearSolver<
      GridType,
      VpeStencilGeneratorType,
      VxpeRhsGeneratorType,
      VxpeResultAcquirerType>
    VxpeSolverType;

  typedef FluidSimulation::LinearSolver<
      GridType,
      VpeStencilGeneratorType,
      VypeRhsGeneratorType,
      VypeResultAcquirerType>
    VypeSolverType;

  typedef typename GhostLayer::Handlers<TD> GhostHandlersType;

  typedef
    VtkPlot<TMemory, TGridGeometry, TScalar, TD>
    VtkPlotType;

  using Plotter = XdmfHdf5Output::XdmfHdf5Writer<GridType>;

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

    _vpeStencilGenerator.initialize(&_parallelDistribution,
                                    &_parameters,
                                    &_dt);
    _ppeStencilGenerator.initialize(&_parallelDistribution);
    _ppeRhsGenerator.initialize(&_dt);
    _ppeResultAcquierer2.initialize(&_dt);

    _vxpeSolver.initialize(&_grid,
                           &_parallelDistribution,
                           &_vpeStencilGenerator,
                           &_vxpeRhsGenerator,
                           &_vxpeResultGenerator,
                           &_ghostHandler.vxpeStencilGeneratorStack,
                           &_ghostHandler.vxpeRhsGeneratorStack,
                           &_ghostHandler.vxpeRhsAcquiererStack);

    _vypeSolver.initialize(&_grid,
                           &_parallelDistribution,
                           &_vpeStencilGenerator,
                           &_vypeRhsGenerator,
                           &_vypeResultGenerator,
                           &_ghostHandler.vypeStencilGeneratorStack,
                           &_ghostHandler.vypeRhsGeneratorStack,
                           &_ghostHandler.vypeRhsAcquiererStack);

    _ppeSolver1.initialize(&_grid,
                           &_parallelDistribution,
                           &_ppeStencilGenerator,
                           &_ppeRhsGenerator,
                           &_ppeResultAcquierer1,
                           &_ghostHandler.ppeStencilGeneratorStack,
                           &_ghostHandler.ppeRhsGeneratorStack,
                           &_ghostHandler.ppeRhsAcquiererStack);

    _ppeSolver2.initialize(&_grid,
                           &_parallelDistribution,
                           &_ppeStencilGenerator,
                           &_ppeRhsGenerator,
                           &_ppeResultAcquierer2,
                           &_ghostHandler.ppeStencilGeneratorStack,
                           &_ghostHandler.ppeRhsGeneratorStack,
                           &_ghostHandler.ppeRhsAcquiererStack);

    _plotter.initialize(&_parallelDistribution,
                        &_grid,
                        outputDirectory,
                        fileNamePrefix);

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
      _plotter.writeGeometry();
      _plot.initialize(&_grid,
                       &_parallelDistribution,
                       &_gridGeometry,
                       outputDirectory,
                       fileNamePrefix);

      _plotter.writeAttributes(_iterationCount);
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

      VelocityType temp = VelocityType::Zero();

      VelocityType velocity;

      for (int d = 0; d < TD; ++d) {
        velocity(d) = 0.5 * (accessor.leftCellInDimension(d)->velocity(d)
                             + accessor.currentCell()->velocity(d));
      }

      accessor.currentCell()->convection() = convection;

      if (TSolverType == 1) {
        auto pressureGradient = computeProjectionPressureGradient(accessor);

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
        accessor.currentCell()->fgh()
          = accessor.currentCell()->velocity()
            + _dt * (-1.5 * convection + 0.5 * previousConvection
                     + 0.5 * diffusion
                     - pressureGradient
                     + _parameters.g());
      } else {
        auto pressureGradient = computePressureGradient(accessor);

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
      _vxpeSolver.update();
      _vxpeSolver.solve();
      _vypeSolver.update();
      _vypeSolver.solve();
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

    if (TSolverType == 1) {
      _ppeSolver2.solve();
    } else {
      _ppeSolver1.solve();
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandler.mpiPressureExchangeStack[d][d2](PETSC_COMM_WORLD);
      }
    }

    for (auto const& accessor : _grid.innerGrid) {
      for (int d = 0; d < TD; ++d) {
        if (TSolverType == 1) {
          accessor.currentCell()->velocity() (d)
            = accessor.currentCell()->fgh() (d)
              - _dt / (0.5 * (accessor.rightWidthInDimension(d)(d)
                              + accessor.currentWidth() (d)))
              * (accessor.rightCellInDimension(d)->pressureProjection()
                 - accessor.currentCell()->pressureProjection());
        } else {
          accessor.currentCell()->velocity() (d)
            = accessor.currentCell()->fgh() (d)
              - _dt / (0.5 * (accessor.rightWidthInDimension(d)(d)
                              + accessor.currentWidth() (d)))
              * (accessor.rightCellInDimension(d)->pressure()
                 - accessor.currentCell()->pressure());
        }
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

        _plotter.writeAttributes(_iterationCount);
        //_plot.plot(_iterationCount, _time, _dt);
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
  VpeStencilGeneratorType                     _vpeStencilGenerator;
  VxpeRhsGeneratorType                        _vxpeRhsGenerator;
  VxpeResultAcquirerType                      _vxpeResultGenerator;
  VypeRhsGeneratorType                        _vypeRhsGenerator;
  VypeResultAcquirerType                      _vypeResultGenerator;
  VxpeSolverType                              _vxpeSolver;
  VypeSolverType                              _vypeSolver;
  PpeStencilGeneratorType                     _ppeStencilGenerator;
  PpeRhsGeneratorType                         _ppeRhsGenerator;
  PpeResultAcquirer1Type                      _ppeResultAcquierer1;
  PpeResultAcquirer2Type                      _ppeResultAcquierer2;
  PpeSolver1Type                              _ppeSolver1;
  PpeSolver2Type                              _ppeSolver2;
  VelocityType                                _maxVelocity;
  GhostHandlersType                           _ghostHandler;
  VtkPlotType                                 _plot;
  Plotter                                     _plotter;
};
}
}
#endif
