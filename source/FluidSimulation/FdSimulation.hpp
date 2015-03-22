#ifndef FsiSimulation_FluidSimulation_FdSimulation_hpp
#define FsiSimulation_FluidSimulation_FdSimulation_hpp

#include "Configuration.hpp"
#include "GhostLayer/Handlers.hpp"
#include "Grid.hpp"
#include "ImmersedBoundary/BodyForce/functions.hpp"
#include "ImmersedBoundary/Fadlun/functions.hpp"
#include "ImmersedBoundary/FeedbackForcing/functions.hpp"
#include "LinearSolver.hpp"
#include "Memory.hpp"
#include "Parameters.hpp"
#include "Private/mpigenerics.hpp"
#include "Simulation.hpp"
#include "functions.hpp"

#include <Uni/Logging/macros>
#include <Uni/Stopwatch>

#include <map>
#include <unordered_map>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TGridGeometry,
          typename TSimulationResultsWriter,
          typename TScalar,
          int TD,
          int TSolverType = 0>
class FdSimulation : public Simulation {
public:
  using BaseType = Simulation;

  using Path = typename BaseType::Path;

  using ParametersType = Parameters<TScalar, TD>;

  using GridGeometryType = TGridGeometry;

  using MemoryType = Memory<GridGeometryType, ParametersType>;

  using CellAccessorType
          = typename MemoryType::CellAccessorType;

  using GridType
          = typename MemoryType::GridType;

  using ParallelDistributionType
          = typename MemoryType::ParallelDistributionType;

  using GhostHandlersType
          = typename GhostLayer::Handlers<Dimensions>;


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

public:
  FdSimulation() {}

  FdSimulation(FdSimulation const& other) = delete;

  ~FdSimulation() {}

  FdSimulation const&
  operator=(FdSimulation const& other) = delete;

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
    int position = accessor.positionInRespectToGeometry();

    bool result = true;

    for (int currentDistance = 1;
         currentDistance <= distance;
         ++currentDistance) {
      for (int d = 0; d < TD; ++d) {
        for (int d2 = 0; d2 < 2; ++d2) {
          auto index = accessor.index();

          if (d2 == 0) {
            index(d) -= currentDistance;
          } else {
            index(d) += currentDistance;
          }

          if ((index(d) < _memory.grid()->innerGrid.leftIndent() (d))
              || (index(d) >= _memory.grid()->innerGrid.innerLimit() (d))) {
            continue;
          }

          auto const position1 = accessor.absoluteCell(index)->positionInRespectToGeometry();

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

    _vpeStencilGenerator.initialize(_memory.parallelDistribution(),
                                    _memory.parameters(),
                                    &_memory.timeStepSize());
    _ppeStencilGenerator.initialize(_memory.parallelDistribution());
    _ppeRhsGenerator.initialize(&_memory.timeStepSize());
    _ppeResultAcquierer2.initialize(&_memory.timeStepSize());

    _vxpeSolver.initialize(_memory.grid(),
                           _memory.parallelDistribution(),
                           &_vpeStencilGenerator,
                           &_vxpeRhsGenerator,
                           &_vxpeResultGenerator,
                           &_ghostHandlers.vpeStencilGeneratorStack[0],
                           &_ghostHandlers.vpeRhsGeneratorStack[0],
                           &_ghostHandlers.vpeRhsAcquiererStack[0]);

    _vypeSolver.initialize(_memory.grid(),
                           _memory.parallelDistribution(),
                           &_vpeStencilGenerator,
                           &_vypeRhsGenerator,
                           &_vypeResultGenerator,
                           &_ghostHandlers.vpeStencilGeneratorStack[1],
                           &_ghostHandlers.vpeRhsGeneratorStack[1],
                           &_ghostHandlers.vpeRhsAcquiererStack[1]);

    _ppeSolver1.initialize(_memory.grid(),
                           _memory.parallelDistribution(),
                           &_ppeStencilGenerator,
                           &_ppeRhsGenerator,
                           &_ppeResultAcquierer1,
                           &_ghostHandlers.ppeStencilGeneratorStack,
                           &_ghostHandlers.ppeRhsGeneratorStack,
                           &_ghostHandlers.ppeRhsAcquiererStack);

    _ppeSolver2.initialize(_memory.grid(),
                           _memory.parallelDistribution(),
                           &_ppeStencilGenerator,
                           &_ppeRhsGenerator,
                           &_ppeResultAcquierer2,
                           &_ghostHandlers.ppeStencilGeneratorStack,
                           &_ghostHandlers.ppeRhsGeneratorStack,
                           &_ghostHandlers.ppeRhsAcquiererStack);

    _resultWriter.initialize(_memory.parallelDistribution(),
                             _memory.grid(),
                             outputDirectory,
                             fileNamePrefix);

    _preciceInterface->initialize();
    _preciceInterface->initializeData();

    auto fluidMeshId = _preciceInterface->getMeshID("FluidMesh");

    for (auto const& accessor : _memory.grid()->innerGrid) {
      ImmersedBoundary::Fadlun::template
      computeDistances<CellAccessorType, TD>(accessor, _preciceInterface);
    }

    for (auto const& accessor : _memory.grid()->innerGrid) {
      Eigen::Vector4i pattern({ 1, 1, 1, 1 });

      bool doAdd = false;

      int distance = computeCellDistance(accessor, 2);

      if (distance != -1) {
        for (int i = 0; i < 4; ++i) {
          auto const& position = pattern(i);

          if (!((i == distance)
                && (position == 1))) {
            continue;
          }

          if (i < 2) {
            if (isOutside(accessor.positionInRespectToGeometry()) == 1) {
              doAdd = true;
            }
          } else {
            if (isOutside(accessor.positionInRespectToGeometry()) == 0) {
              doAdd = true;
            }
          }
          break;
        }
      }

      if (doAdd) {
        VectorDsType position = accessor.position();
        position += 0.5 * accessor.width();
        auto vertexId =
          _preciceInterface->setMeshVertex(fluidMeshId,
                                           position.data());

        _vertexIds.insert(std::make_pair(accessor.indexValues(), vertexId));
      }
    }

    for (auto const& accessor : * _memory.grid()) {
      accessor.velocity() = VelocityType::Zero();
      accessor.fgh()      = VelocityType::Zero();
      accessor.pressure() = 0.0;
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandlers.velocityInitialization[d][d2]();
      }
    }

    _ghostHandlers.executeVelocityMpiExchange();

    _lastPlotTimeStamp = 0.0;

    _memory.maxVelocity()     = VectorDsType::Zero();
    _memory.timeStepSize()    = 1.0;
    _memory.time()            = 0.0;
    _memory.iterationNumber() = 0;

    if (_plotInterval >= 0) {
      _resultWriter.writeGeometry();
      _resultWriter.writeAttributes(_memory.iterationNumber(), _memory.time());
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
    _memory.maxVelocity() = VelocityType::Zero();

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

      auto diffusion = DiffusionProcessing<TD>::compute(accessor);

      diffusion = (1.0 / _memory.parameters()->re()) * diffusion;

      VelocityType temp = VelocityType::Zero();

      VelocityType velocity;

      for (int d = 0; d < TD; ++d) {
        velocity(d) = 0.5 * (accessor.leftCellInDimension(d)->velocity(d)
                             + accessor.velocity(d));
      }

      accessor.convection() = convection;

      if (TSolverType == 1) {
        auto pressureGradient = computeProjectionPressureGradient(accessor);

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
            + _memory.timeStepSize() * (-1.5 * convection + 0.5 *
                                        previousConvection
                                        + 0.5 * diffusion
                                        - pressureGradient
                                        + _memory.parameters()->g());
      } else {
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
    }

    _preciceInterface->advance(_memory.timeStepSize());

    for (auto const& accessor : _memory.grid()->innerGrid) {
      auto find_it =  _vertexIds.find(accessor.index());

      if (find_it == _vertexIds.end()) {
        continue;
      }
      VelocityType force;
      _preciceInterface->readVectorData(
        fluidMeshForcesId,
        find_it->second,
        force.data());

      accessor.fgh() += force;
    }

    if (TSolverType == 1) {
      _vxpeSolver.update();
      _vxpeSolver.solve();
      _vypeSolver.update();
      _vypeSolver.solve();
    } else {
      for (int d = 0; d < TD; ++d) {
        for (int d2 = 0; d2 < 2; ++d2) {
          _ghostHandlers.fghInitialization[d][d2]();
        }
      }
    }

    _ghostHandlers.executeFghMpiExchange();

    if (TSolverType == 1) {
      _ppeSolver2.solve();
    } else {
      _ppeSolver1.solve();
    }

    _ghostHandlers.executePressureMpiExchange();

    for (auto const& accessor : _memory.grid()->innerGrid) {
      for (int d = 0; d < TD; ++d) {
        if (TSolverType == 1) {
          accessor.velocity() (d)
            = accessor.fgh() (d)
              - _memory.timeStepSize() /
              (0.5 * (accessor.rightWidthInDimension(d)(d)
                      + accessor.width() (d)))
              * (accessor.rightCellInDimension(d)->pressureProjection()
                 - accessor.pressureProjection());
        } else {
          accessor.velocity() (d)
            = accessor.fgh() (d)
              - _memory.timeStepSize() /
              (0.5 * (accessor.rightWidthInDimension(d)(d)
                      + accessor.width() (d)))
              * (accessor.rightCellInDimension(d)->pressure()
                 - accessor.pressure());
        }
      }
      computeMaxVelocity<CellAccessorType, TScalar, TD>
        (accessor, _memory.maxVelocity());
    }

    for (int d = 0; d < TD; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostHandlers.velocityInitialization[d][d2]();
      }
    }

    _ghostHandlers.executeVelocityMpiExchange();

    VectorDsType force = VectorDsType::Zero();

    for (auto const& accessor : _memory.grid()->innerGrid) {
      ImmersedBoundary::BodyForce::
      template computeCellForce<CellAccessorType, TD>(
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

        _resultWriter.writeAttributes(_memory.iterationNumber(),
                                      _memory.time());
      }
    }

    force *= 2 / (0.3 * 0.3 * 0.1);

    if (_parallelDistribution->rank == 0) {
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
  std::unordered_map<VectorDiType, int, Hash> _vertexIds;
  precice::SolverInterface*                   _preciceInterface;
  long double                                 _timeLimit;
  TScalar                                     _lastPlotTimeStamp;
  TScalar                                     _plotInterval;
  unsigned long long                          _iterationLimit;
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
  GhostHandlersType                           _ghostHandlers;
  TSimulationResultsWriter                    _resultWriter;
};
}
}
#endif
