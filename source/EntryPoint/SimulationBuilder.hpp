#ifndef FsiSimulation_EntryPoint_SimulationBuilder_hpp
#define FsiSimulation_EntryPoint_SimulationBuilder_hpp

#include "FluidSimulation/Cell.hpp"
#include "FluidSimulation/Configuration.hpp"
#include "FluidSimulation/FdSimulation.hpp"
#include "FluidSimulation/GhostLayer/InitializationActions.hpp"
#include "FluidSimulation/GhostLayer/PetscExchangeActions.hpp"
#include "FluidSimulation/GridGeometry.hpp"
#include "FluidSimulation/VtkOutput/VtkWriter.hpp"
#include "FluidSimulation/XdmfHdf5Output/XdmfHdf5Writer.hpp"

#include "StructuredMemory/Accessor.hpp"
#include "StructuredMemory/Memory.hpp"

#include <Uni/StructuredGrid/Basic/MultiIndex>
#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace EntryPoint {
template <typename Scalar, int TDimensions, int TSolverType = 0>
class SimulationBuilder {
  static_assert((TDimensions > 1) && (TDimensions < 4),
                "Only 2D and 3D simulations supported");

public:
  enum {
    Dimensions = TDimensions
  };
  typedef FluidSimulation::FdSimulation <
    FluidSimulation::UniformGridGeometry<Scalar, TDimensions>,
    StructuredMemory::IterableMemory<FluidSimulation::Cell<Scalar, TDimensions>,
                                     TDimensions>,
    FluidSimulation::XdmfHdf5Output::XdmfHdf5Writer < FluidSimulation::Grid<
      StructuredMemory::IterableMemory<FluidSimulation::Cell<Scalar,
                                                             TDimensions>,
                                       TDimensions>,
      FluidSimulation::UniformGridGeometry<Scalar, TDimensions>,
      TDimensions >>,
      Scalar,
      TDimensions,
      TSolverType> Simulation;
  typedef typename Simulation::GridType                 GridS;
  typedef typename Simulation::GridGeometryType         GridGeometryS;
  typedef typename Simulation::ParallelDistributionType ParallelTopologyS;
  typedef typename Simulation::GhostHandlersType        GhostCellsHandlerS;
  typedef typename Simulation::MemoryType               MemoryS;
  typedef typename Simulation::CellAccessorFactory
    CellAccessorFactoryS;
  typedef typename GridS::CellAccessorType CellAccessorS;
  typedef typename CellAccessorS::CellType CellS;
  typedef typename CellS::Velocity         VelocityS;
  typedef typename CellS::Pressure         PressureS;

  template <int TDimension, int TDirection>
  struct GhostHandlers {
    typedef
      FluidSimulation::GhostLayer::MpiExchange::Handler
      <Scalar, typename GridS::Base,
       SimulationBuilder::template fghAccessor<TDimension>,
       Scalar, TDimensions, TDimension, TDirection>
      FghMpiExchange;
    typedef
      FluidSimulation::GhostLayer::MpiExchange::Handler
      <PressureS, typename GridS::Base,
       SimulationBuilder::pressureAccessor, Scalar, TDimensions,
       TDimension, TDirection>
      PressureMpiExchange;
    typedef
      FluidSimulation::GhostLayer::MpiExchange::Handler
      <VelocityS, typename GridS::Base,
       SimulationBuilder::velocityAccessor, Scalar, TDimensions,
       TDimension, TDirection>
      VelocityMpiExchange;

    typedef
      FluidSimulation::GhostLayer::Initialization::MovingWallFghAction
      <typename GridS::Base, Scalar, TDimensions, TDimension, TDirection>
      MovingWallFghInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, MovingWallFghInitializationAction, TDimensions,
       TDimension, TDirection>
      MovingWallFghInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::InputFghAction
      <typename GridS::Base, Scalar, TDimensions, TDimension, TDirection>
      InputFghInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, InputFghInitializationAction, TDimensions,
       TDimension, TDirection>
      InputFghInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::ParabolicInputFghAction
      <typename GridS::Base, Scalar, TDimensions, TDimension, TDirection>
      ParabolicInputFghInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, ParabolicInputFghInitializationAction, TDimensions,
       TDimension, TDirection>
      ParabolicInputFghInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::OutputFghAction
      <typename GridS::Base, Scalar, TDimensions, TDimension, TDirection>
      OutputFghInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, OutputFghInitializationAction, TDimensions,
       TDimension, TDirection>
      OutputFghInitialization;

    typedef
      FluidSimulation::GhostLayer::LsStencilGenerator::Handler
      <GridS, TDimension, TDirection>
      PpeStencilGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::ConstantRhsGenerationAction<
        CellAccessorS, TDimension>
      PpeRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, PpeRhsGenerationAction, TDimension, TDirection>
      PpeRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::PpeRhsAcquiererAction
      PpeRhsAcquiererAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, PpeRhsAcquiererAction, TDimension, TDirection>
      PpeRhsAcquiererHandler;

    typedef
      FluidSimulation::GhostLayer::LsStencilGenerator::Handler
      <GridS, TDimension, TDirection>
      VpeStencilGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::ConstantRhsGenerationAction
      <CellAccessorS, TDimension>
      VpeConstantRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VpeConstantRhsGenerationAction, TDimension, TDirection>
      VpeConstantRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeInputRhsGenerationAction
      <0, TDimension, TDirection>
      VxpeInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VxpeInputRhsGenerationAction, TDimension, TDirection>
      VxpeInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::
      VpeParabolicInputRhsGenerationAction
      <0, TDimension, TDirection>
      VxpeParabolicInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VxpeParabolicInputRhsGenerationAction, TDimension, TDirection>
      VxpeParabolicInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeRhsAcquiererAction<0>
      VxpeRhsAcquiererAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VxpeRhsAcquiererAction, TDimension, TDirection>
      VxpeRhsAcquiererHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeInputRhsGenerationAction
      <1, TDimension, TDirection>
      VypeInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VypeInputRhsGenerationAction, TDimension, TDirection>
      VypeInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::
      VpeParabolicInputRhsGenerationAction
      <1, TDimension, TDirection>
      VypeParabolicInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VypeParabolicInputRhsGenerationAction, TDimension, TDirection>
      VypeParabolicInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeRhsAcquiererAction<0>
      VypeRhsAcquiererAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VypeRhsAcquiererAction, TDimension, TDirection>
      VypeRhsAcquiererHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeInputRhsGenerationAction
      <2, TDimension, TDirection>
      VzpeInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VzpeInputRhsGenerationAction, TDimension, TDirection>
      VzpeInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::
      VpeParabolicInputRhsGenerationAction
      <2, TDimension, TDirection>
      VzpeParabolicInputRhsGenerationAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VzpeParabolicInputRhsGenerationAction, TDimension, TDirection>
      VzpeParabolicInputRhsGenerationHandler;

    typedef
      FluidSimulation::GhostLayer::PetscExchange::VpeRhsAcquiererAction<0>
      VzpeRhsAcquiererAction;
    typedef
      FluidSimulation::GhostLayer::PetscExchange::Handler
      <GridS, VzpeRhsAcquiererAction, TDimension, TDirection>
      VzpeRhsAcquiererHandler;

    typedef
      FluidSimulation::GhostLayer::Initialization::MovingWallVelocityAction
      <typename GridS::Base, Scalar, TDimensions, TDimension, TDirection>
      MovingWallVelocityInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, MovingWallVelocityInitializationAction,
       TDimensions,
       TDimension, TDirection>
      MovingWallVelocityInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::InputVelocityAction
      <typename GridS::Base, Scalar, TDimensions, TDimension, TDirection>
      InputVelocityInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, InputVelocityInitializationAction, TDimensions,
       TDimension, TDirection>
      InputVelocityInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::ParabolicInputVelocityAction
      <typename GridS::Base, Scalar, TDimensions, TDimension, TDirection>
      ParabolicInputVelocityInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, ParabolicInputVelocityInitializationAction,
       TDimensions,
       TDimension, TDirection>
      ParabolicInputVelocityInitialization;
    typedef
      FluidSimulation::GhostLayer::Initialization::OutputVelocityAction
      <typename GridS::Base, Scalar, TDimensions, TDimension, TDirection>
      OutputVelocityInitializationAction;
    typedef
      FluidSimulation::GhostLayer::Initialization::Handler
      <typename GridS::Base, OutputVelocityInitializationAction, TDimensions,
       TDimension, TDirection>
      OutputVelocityInitialization;

public:
    GhostHandlers(
      FluidSimulation::Configuration* configuration,
      GridS*                          grid,
      ParallelTopologyS*              parallelDistribution,
      GhostCellsHandlerS*             handlers,
      VelocityS*                      maxVelocity)
      : _configuration(configuration),
      _grid(grid),
      _parallelDistribution(parallelDistribution),
      _handlers(handlers),
      _maxVelocity(maxVelocity) {}

    inline void
    setAsInput() {
      _handlers->fghInitialization[TDimension][TDirection] =
        InputFghInitialization::getHandler(
          &_grid->indentedBoundaries[TDimension][TDirection],
          _parallelDistribution,
          new InputFghInitializationAction(_configuration));

      _handlers->ppeStencilGeneratorStack[TDimension][TDirection] =
        PpeStencilGenerationHandler::getNeumannMiddle(_grid,
                                                      _parallelDistribution);
      _handlers->ppeRhsGeneratorStack[TDimension][TDirection] =
        PpeRhsGenerationHandler::getHandler(_grid, _parallelDistribution,
                                            new PpeRhsGenerationAction(0.0));

      if (TDirection == 1) {
        _handlers->ppeRhsAcquiererStack[TDimension][TDirection] =
          PpeRhsAcquiererHandler::getHandler(
            _grid, _parallelDistribution, new PpeRhsAcquiererAction());
      }

      for (int d = 0; d < Dimensions; ++d) {
        if (d == TDimension) {
          if (TDirection == 0) {
            _handlers->vpeStencilGeneratorStack[d][TDimension][TDirection] =
              VpeStencilGenerationHandler::getDirichletLeft(_grid,
                                                            _parallelDistribution);
          } else {
            _handlers->vpeStencilGeneratorStack[d][TDimension][TDirection] =
              VpeStencilGenerationHandler::getDirichletRight(_grid,
                                                             _parallelDistribution);
          }
        } else {
          _handlers->vpeStencilGeneratorStack[d][TDimension][TDirection] =
            VpeStencilGenerationHandler::getDirichletMiddle(_grid,
                                                            _parallelDistribution);
        }
      }

      _handlers->vpeRhsGeneratorStack[0][TDimension][TDirection] =
        VxpeInputRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VxpeInputRhsGenerationAction(_configuration));
      _handlers->vpeRhsAcquiererStack[0][TDimension][TDirection] =
        VxpeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VxpeRhsAcquiererAction());

      _handlers->vpeRhsGeneratorStack[1][TDimension][TDirection] =
        VypeInputRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VypeInputRhsGenerationAction(_configuration));
      _handlers->vpeRhsAcquiererStack[1][TDimension][TDirection] =
        VypeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VypeRhsAcquiererAction());

      if (Dimensions > 2) {
        _handlers->vpeRhsGeneratorStack[2][TDimension][TDirection] =
          VzpeInputRhsGenerationHandler::getHandler(
            _grid,
            _parallelDistribution,
            new VzpeInputRhsGenerationAction(_configuration));
        _handlers->vpeRhsAcquiererStack[2][TDimension][TDirection] =
          VzpeRhsAcquiererHandler::getHandler(
            _grid, _parallelDistribution, new VzpeRhsAcquiererAction());
      }

      _handlers->velocityInitialization[TDimension][TDirection] =
        InputVelocityInitialization::getHandler(
          &_grid->boundaries[TDimension][TDirection],
          _parallelDistribution,
          new InputVelocityInitializationAction(_configuration, _maxVelocity));
    }

    inline void
    setAsParabolicInput() {
      _handlers->fghInitialization[TDimension][TDirection] =
        ParabolicInputFghInitialization::getHandler(
          &_grid->indentedBoundaries[TDimension][TDirection],
          _parallelDistribution,
          new ParabolicInputFghInitializationAction(_configuration));

      _handlers->ppeStencilGeneratorStack[TDimension][TDirection] =
        PpeStencilGenerationHandler::getNeumannMiddle(_grid,
                                                      _parallelDistribution);
      _handlers->ppeRhsGeneratorStack[TDimension][TDirection] =
        PpeRhsGenerationHandler::getHandler(_grid, _parallelDistribution,
                                            new PpeRhsGenerationAction(0.0));

      if (TDirection == 1) {
        _handlers->ppeRhsAcquiererStack[TDimension][TDirection] =
          PpeRhsAcquiererHandler::getHandler(
            _grid, _parallelDistribution, new PpeRhsAcquiererAction());
      }

      for (int d = 0; d < Dimensions; ++d) {
        if (d == TDimension) {
          if (TDirection == 0) {
            _handlers->vpeStencilGeneratorStack[d][TDimension][TDirection] =
              VpeStencilGenerationHandler::getDirichletLeft(_grid,
                                                            _parallelDistribution);
          } else {
            _handlers->vpeStencilGeneratorStack[d][TDimension][TDirection] =
              VpeStencilGenerationHandler::getDirichletRight(_grid,
                                                             _parallelDistribution);
          }
        } else {
          _handlers->vpeStencilGeneratorStack[d][TDimension][TDirection] =
            VpeStencilGenerationHandler::getDirichletMiddle(_grid,
                                                            _parallelDistribution);
        }
      }

      _handlers->vpeRhsGeneratorStack[0][TDimension][TDirection] =
        VxpeParabolicInputRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VxpeParabolicInputRhsGenerationAction(_configuration));
      _handlers->vpeRhsAcquiererStack[0][TDimension][TDirection] =
        VxpeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VxpeRhsAcquiererAction());

      _handlers->vpeRhsGeneratorStack[1][TDimension][TDirection] =
        VypeParabolicInputRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VypeParabolicInputRhsGenerationAction(_configuration));
      _handlers->vpeRhsAcquiererStack[1][TDimension][TDirection] =
        VypeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VypeRhsAcquiererAction());

      if (Dimensions > 2) {
        _handlers->vpeRhsGeneratorStack[2][TDimension][TDirection] =
          VzpeInputRhsGenerationHandler::getHandler(
            _grid,
            _parallelDistribution,
            new VzpeInputRhsGenerationAction(_configuration));
        _handlers->vpeRhsAcquiererStack[2][TDimension][TDirection] =
          VzpeRhsAcquiererHandler::getHandler(
            _grid, _parallelDistribution, new VzpeRhsAcquiererAction());
      }

      _handlers->velocityInitialization[TDimension][TDirection] =
        ParabolicInputVelocityInitialization::getHandler(
          &_grid->boundaries[TDimension][TDirection],
          _parallelDistribution,
          new ParabolicInputVelocityInitializationAction(_configuration,
                                                         _maxVelocity));
    }

    inline void
    setAsOutput() {
      _handlers->fghInitialization[TDimension][TDirection] =
        OutputFghInitialization::getHandler(
          &_grid->indentedBoundaries[TDimension][TDirection],
          _parallelDistribution,
          new OutputFghInitializationAction());

      _handlers->ppeStencilGeneratorStack[TDimension][TDirection] =
        PpeStencilGenerationHandler::getDirichletMiddle(_grid,
                                                        _parallelDistribution);
      _handlers->ppeRhsGeneratorStack[TDimension][TDirection] =
        PpeRhsGenerationHandler::getHandler(_grid, _parallelDistribution,
                                            new PpeRhsGenerationAction(0.0));

      if (TDirection == 1) {
        _handlers->ppeRhsAcquiererStack[TDimension][TDirection] =
          PpeRhsAcquiererHandler::getHandler(
            _grid, _parallelDistribution, new PpeRhsAcquiererAction());
      }

      for (int d = 0; d < Dimensions; ++d) {
        if (d == TDimension) {
          if (TDirection == 0) {
            _handlers->vpeStencilGeneratorStack[d][TDimension][TDirection] =
              VpeStencilGenerationHandler::getNeumannLeft(_grid,
                                                          _parallelDistribution);
          } else {
            _handlers->vpeStencilGeneratorStack[d][TDimension][TDirection] =
              VpeStencilGenerationHandler::getNeumannRight(_grid,
                                                           _parallelDistribution);
          }
        } else {
          _handlers->vpeStencilGeneratorStack[d][TDimension][TDirection] =
            VpeStencilGenerationHandler::getNeumannMiddle(_grid,
                                                          _parallelDistribution);
        }
      }

      int offset = 0;

      if (TDimension == 0) {
        if (TDirection == 1) {
          offset = 1;
        }
      }
      _handlers->vpeRhsGeneratorStack[0][TDimension][TDirection] =
        VpeConstantRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VpeConstantRhsGenerationAction(0.0, offset));
      _handlers->vpeRhsAcquiererStack[0][TDimension][TDirection] =
        VxpeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VxpeRhsAcquiererAction());

      offset = 0;

      if (TDimension == 1) {
        if (TDirection == 1) {
          offset = 1;
        }
      }

      _handlers->vpeRhsGeneratorStack[1][TDimension][TDirection] =
        VpeConstantRhsGenerationHandler::getHandler(
          _grid,
          _parallelDistribution,
          new VpeConstantRhsGenerationAction(0.0, offset));
      _handlers->vpeRhsAcquiererStack[1][TDimension][TDirection] =
        VypeRhsAcquiererHandler::getHandler(
          _grid, _parallelDistribution, new VypeRhsAcquiererAction());

      if (Dimensions > 2) {
        offset = 0;

        if (TDimension == 2) {
          if (TDirection == 1) {
            offset = 1;
          }
        }
        _handlers->vpeRhsGeneratorStack[2][TDimension][TDirection] =
          VpeConstantRhsGenerationHandler::getHandler(
            _grid,
            _parallelDistribution,
            new VpeConstantRhsGenerationAction(0.0, offset));
        _handlers->vpeRhsAcquiererStack[2][TDimension][TDirection] =
          VzpeRhsAcquiererHandler::getHandler(
            _grid, _parallelDistribution, new VzpeRhsAcquiererAction());
      }

      _handlers->velocityInitialization[TDimension][TDirection] =
        OutputVelocityInitialization::getHandler(
          &_grid->boundaries[TDimension][TDirection],
          _parallelDistribution,
          new OutputVelocityInitializationAction(_maxVelocity));
    }

private:
    FluidSimulation::Configuration* _configuration;
    GridS*                          _grid;
    ParallelTopologyS*              _parallelDistribution;
    GhostCellsHandlerS*             _handlers;
    VelocityS*                      _maxVelocity;
  };

public:
  SimulationBuilder(FluidSimulation::Configuration* configuration)
    : _configuration(configuration) {
    _simulation           = new Simulation();
    _grid                 = &_simulation->_grid;
    _parallelDistribution = &_simulation->_parallelDistribution;
    _handlers             = &_simulation->_ghostHandler;
    _maxVelocity          = &_simulation->_maxVelocity;

    typedef typename ParallelTopologyS::VectorDi VectorDi;
    typedef typename GridGeometryS::VectorDs     VectorDs;

    VectorDi processorSize;
    VectorDi globalCellSize;
    VectorDs width;
    processorSize(0)  = configuration->parallelizationSize(0);
    globalCellSize(0) = configuration->size(0);
    width(0)          = (Scalar)configuration->width(0);
    processorSize(1)  = configuration->parallelizationSize(1);
    globalCellSize(1) = configuration->size(1);
    width(1)          = (Scalar)configuration->width(1);

    if (TDimensions == 3) {
      processorSize(2)  = configuration->parallelizationSize(2);
      globalCellSize(2) = configuration->size(2);
      width(2)          = (Scalar)configuration->width(2);
    }

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    _parallelDistribution->initialize(rank,
                                      processorSize,
                                      globalCellSize);

    _simulation->_parameters.re()    = configuration->re;
    _simulation->_parameters.gamma() = configuration->gamma;
    _simulation->_parameters.tau()   = configuration->tau;
    _simulation->_parameters.alpha() = configuration->alpha;
    _simulation->_parameters.immersedBoundaryMethod()
                                  = configuration->immersedBoundaryMethod;
    _simulation->_parameters.g(0) = configuration->environment(0);
    _simulation->_parameters.g(1) = configuration->environment(1);
    _simulation->_iterationLimit  = configuration->iterationLimit;
    _simulation->_timeLimit       = configuration->timeLimit;
    _simulation->_plotInterval    = configuration->plotInterval;

    if (TDimensions == 3) {
      _simulation->_parameters.g(2) = configuration->environment(2);
    }

    VectorDi localSize(_parallelDistribution->localCellSize + 2 *
                       VectorDi::Ones());

    _simulation->setParameters(localSize, width);

    setLeftWallAs(configuration->walls[0][0]->type());
    setRightWallAs(configuration->walls[0][1]->type());
    setBottomWallAs(configuration->walls[1][0]->type());
    setTopWallAs(configuration->walls[1][1]->type());
    setBackWallAs(configuration->walls[2][0]->type());
    setFrontWallAs(configuration->walls[2][1]->type());
  }

  Simulation*
  simulation() const {
    return _simulation;
  }

  void
  setLeftWallAs(FluidSimulation::WallType type) {
    if (_parallelDistribution->neighbors[0][0] >= 0) {
      _setLeftAsMpiExchange();
    } else if (type == FluidSimulation::WallType::Input) {
      _setLeftAsInput();
    } else if (type == FluidSimulation::WallType::ParabolicInput) {
      _setLeftAsParabolicInput();
    } else if (type == FluidSimulation::WallType::Output) {
      _setLeftAsOutput();
    }
  }

  void
  setRightWallAs(FluidSimulation::WallType type) {
    if (_parallelDistribution->neighbors[0][1] >= 0) {
      _setRightAsMpiExchange();
    } else if (type == FluidSimulation::WallType::Input) {
      _setRightAsInput();
    } else if (type == FluidSimulation::WallType::ParabolicInput) {
      _setRightAsParabolicInput();
    } else if (type == FluidSimulation::WallType::Output) {
      _setRightAsOutput();
    }
  }

  void
  setBottomWallAs(FluidSimulation::WallType type) {
    if (Dimensions < 2) { return; }

    if (_parallelDistribution->neighbors[1][0] >= 0) {
      _setBottomAsMpiExchange();
    } else if (type == FluidSimulation::WallType::Input) {
      _setBottomAsInput();
    } else if (type == FluidSimulation::WallType::ParabolicInput) {
      _setBottomAsParabolicInput();
    } else if (type == FluidSimulation::WallType::Output) {
      _setBottomAsOutput();
    }
  }

  void
  setTopWallAs(FluidSimulation::WallType type) {
    if (Dimensions < 2) { return; }

    if (_parallelDistribution->neighbors[1][1] >= 0) {
      _setTopAsMpiExchange();
    } else if (type == FluidSimulation::WallType::Input) {
      _setTopAsInput();
    } else if (type == FluidSimulation::WallType::ParabolicInput) {
      _setTopAsParabolicInput();
    } else if (type == FluidSimulation::WallType::Output) {
      _setTopAsOutput();
    }
  }

  void
  setBackWallAs(FluidSimulation::WallType type) {
    if (Dimensions < 3) { return; }

    if (_parallelDistribution->neighbors[2][0] >= 0) {
      _setBackAsMpiExchange();
    } else if (type == FluidSimulation::WallType::Input) {
      _setBackAsInput();
    } else if (type == FluidSimulation::WallType::ParabolicInput) {
      _setBackAsParabolicInput();
    } else if (type == FluidSimulation::WallType::Output) {
      _setBackAsOutput();
    }
  }

  void
  setFrontWallAs(FluidSimulation::WallType type) {
    if (Dimensions < 3) { return; }

    if (_parallelDistribution->neighbors[2][1] >= 0) {
      _setFrontAsMpiExchange();
    } else if (type == FluidSimulation::WallType::Input) {
      _setFrontAsInput();
    } else if (type == FluidSimulation::WallType::ParabolicInput) {
      _setFrontAsParabolicInput();
    } else if (type == FluidSimulation::WallType::Output) {
      _setFrontAsOutput();
    }
  }

private:
  void
  _setLeftAsInput() {
    typedef GhostHandlers<0, 0> LeftHandlers;
    typedef typename LeftHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    LeftHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsInput();
  }

  void
  _setLeftAsParabolicInput() {
    typedef GhostHandlers<0, 0> LeftHandlers;
    typedef typename LeftHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    LeftHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsParabolicInput();
  }

  void
  _setLeftAsOutput() {
    typedef GhostHandlers<0, 0> LeftHandlers;

    typedef typename LeftHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    LeftHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsOutput();
  }

  void
  _setLeftAsMpiExchange() {
    typedef GhostHandlers<0, 0>                        LeftHandlers;
    typedef typename LeftHandlers::FghMpiExchange      FghHandler;
    typedef typename LeftHandlers::PressureMpiExchange PressureHandler;
    typedef typename LeftHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[0][0] =
      FghHandler::getReceiveHandler(_grid,
                                    _parallelDistribution,
                                    leftIndent,
                                    rightIndent);

    _handlers->mpiPressureExchangeStack[0][0] =
      PressureHandler::getSendHandler(_grid,
                                      _parallelDistribution,
                                      leftIndent,
                                      rightIndent);

    rightIndent(1) = 0;
    rightIndent(2) = 0;

    _handlers->mpiVelocityExchangeStack[0][0] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setRightAsInput() {
    typedef GhostHandlers<0, 1> RightHandlers;

    typedef typename RightHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    RightHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsInput();
  }

  void
  _setRightAsParabolicInput() {
    typedef GhostHandlers<0, 1> RightHandlers;
    typedef typename RightHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    RightHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsParabolicInput();
  }

  void
  _setRightAsOutput() {
    typedef GhostHandlers<0, 1> RightHandlers;

    typedef typename RightHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    RightHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsOutput();
  }

  void
  _setRightAsMpiExchange() {
    typedef GhostHandlers<0, 1>                         RightHandlers;
    typedef typename RightHandlers::FghMpiExchange      FghHandler;
    typedef typename RightHandlers::PressureMpiExchange PressureHandler;
    typedef typename RightHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[0][1] =
      FghHandler::getSendHandler(_grid,
                                 _parallelDistribution,
                                 leftIndent,
                                 rightIndent);

    _handlers->mpiPressureExchangeStack[0][1] =
      PressureHandler::getReceiveHandler(_grid,
                                         _parallelDistribution,
                                         leftIndent,
                                         rightIndent);

    leftIndent(1) = 0;
    leftIndent(2) = 0;

    _handlers->mpiVelocityExchangeStack[0][1] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setBottomAsInput() {
    typedef GhostHandlers<1, 0> BottomHandlers;
    typedef typename BottomHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    BottomHandlers handlers(_configuration,
                            _grid,
                            _parallelDistribution,
                            _handlers,
                            _maxVelocity);

    handlers.setAsInput();
  }

  void
  _setBottomAsParabolicInput() {
    typedef GhostHandlers<1, 0> BottomHandlers;
    typedef typename BottomHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    BottomHandlers handlers(_configuration,
                            _grid,
                            _parallelDistribution,
                            _handlers,
                            _maxVelocity);

    handlers.setAsParabolicInput();
  }

  void
  _setBottomAsOutput() {
    typedef GhostHandlers<1, 0> BottomHandlers;

    typedef typename BottomHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    BottomHandlers handlers(_configuration,
                            _grid,
                            _parallelDistribution,
                            _handlers,
                            _maxVelocity);

    handlers.setAsOutput();
  }

  void
  _setBottomAsMpiExchange() {
    typedef GhostHandlers<1, 0>                          BottomHandlers;
    typedef typename BottomHandlers::FghMpiExchange      FghHandler;
    typedef typename BottomHandlers::PressureMpiExchange PressureHandler;
    typedef typename BottomHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[1][0] =
      FghHandler::getReceiveHandler(_grid,
                                    _parallelDistribution,
                                    leftIndent,
                                    rightIndent);

    _handlers->mpiPressureExchangeStack[1][0] =
      PressureHandler::getSendHandler(_grid,
                                      _parallelDistribution,
                                      leftIndent,
                                      rightIndent);

    rightIndent(0) = 0;
    rightIndent(2) = 0;

    _handlers->mpiVelocityExchangeStack[1][0] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setTopAsInput() {
    typedef GhostHandlers<1, 1> TopHandlers;

    typedef typename TopHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    TopHandlers handlers(_configuration,
                         _grid,
                         _parallelDistribution,
                         _handlers,
                         _maxVelocity);

    handlers.setAsInput();
  }

  void
  _setTopAsParabolicInput() {
    typedef GhostHandlers<1, 1> TopHandlers;
    typedef typename TopHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    TopHandlers handlers(_configuration,
                         _grid,
                         _parallelDistribution,
                         _handlers,
                         _maxVelocity);

    handlers.setAsParabolicInput();
  }

  void
  _setTopAsOutput() {
    typedef GhostHandlers<1, 1> TopHandlers;

    typedef typename TopHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    TopHandlers handlers(_configuration,
                         _grid,
                         _parallelDistribution,
                         _handlers,
                         _maxVelocity);

    handlers.setAsOutput();
  }

  void
  _setTopAsMpiExchange() {
    typedef GhostHandlers<1, 1>                       TopHandlers;
    typedef typename TopHandlers::FghMpiExchange      FghHandler;
    typedef typename TopHandlers::PressureMpiExchange PressureHandler;
    typedef typename TopHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[1][1] =
      FghHandler::getSendHandler(_grid,
                                 _parallelDistribution,
                                 leftIndent,
                                 rightIndent);

    _handlers->mpiPressureExchangeStack[1][1] =
      PressureHandler::getReceiveHandler(_grid,
                                         _parallelDistribution,
                                         leftIndent,
                                         rightIndent);

    leftIndent(0) = 0;
    leftIndent(2) = 0;

    _handlers->mpiVelocityExchangeStack[1][1] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setBackAsInput() {
    typedef GhostHandlers<2, 0> BackHandlers;
    typedef typename BackHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    BackHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsInput();
  }

  void
  _setBackAsParabolicInput() {
    typedef GhostHandlers<2, 0> BackHandlers;
    typedef typename BackHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    BackHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsParabolicInput();
  }

  void
  _setBackAsOutput() {
    typedef GhostHandlers<2, 0> BackHandlers;

    typedef typename BackHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    BackHandlers handlers(_configuration,
                          _grid,
                          _parallelDistribution,
                          _handlers,
                          _maxVelocity);

    handlers.setAsOutput();
  }

  void
  _setBackAsMpiExchange() {
    typedef GhostHandlers<2, 0>                        BackHandlers;
    typedef typename BackHandlers::FghMpiExchange      FghHandler;
    typedef typename BackHandlers::PressureMpiExchange PressureHandler;
    typedef typename BackHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[2][0] =
      FghHandler::getReceiveHandler(_grid,
                                    _parallelDistribution,
                                    leftIndent,
                                    rightIndent);

    _handlers->mpiPressureExchangeStack[2][0] =
      PressureHandler::getSendHandler(_grid,
                                      _parallelDistribution,
                                      leftIndent,
                                      rightIndent);

    rightIndent(0) = 0;
    rightIndent(1) = 0;

    _handlers->mpiVelocityExchangeStack[2][0] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  void
  _setFrontAsInput() {
    typedef GhostHandlers<2, 1> FrontHandlers;

    typedef typename FrontHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    FrontHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsInput();
  }

  void
  _setFrontAsParabolicInput() {
    typedef GhostHandlers<2, 1> FrontHandlers;
    typedef typename FrontHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    FrontHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsParabolicInput();
  }

  void
  _setFrontAsOutput() {
    typedef GhostHandlers<2, 1> FrontHandlers;

    typedef typename FrontHandlers::VpeStencilGenerationHandler
      VpeStencilGenerationHandler;

    FrontHandlers handlers(_configuration,
                           _grid,
                           _parallelDistribution,
                           _handlers,
                           _maxVelocity);

    handlers.setAsOutput();
  }

  void
  _setFrontAsMpiExchange() {
    typedef GhostHandlers<2, 1>                         FrontHandlers;
    typedef typename FrontHandlers::FghMpiExchange      FghHandler;
    typedef typename FrontHandlers::PressureMpiExchange PressureHandler;
    typedef typename FrontHandlers::VelocityMpiExchange VelocityHandler;

    auto leftIndent  = _grid->innerGrid.leftIndent();
    auto rightIndent = _grid->innerGrid.rightIndent();

    _handlers->mpiFghExchangeStack[2][1] =
      FghHandler::getSendHandler(_grid,
                                 _parallelDistribution,
                                 leftIndent,
                                 rightIndent);

    _handlers->mpiPressureExchangeStack[2][1] =
      PressureHandler::getReceiveHandler(_grid,
                                         _parallelDistribution,
                                         leftIndent,
                                         rightIndent);

    leftIndent(0) = 0;
    leftIndent(1) = 0;

    _handlers->mpiVelocityExchangeStack[2][1] =
      VelocityHandler::getExchangeHandler(_grid,
                                          _parallelDistribution,
                                          leftIndent,
                                          rightIndent);
  }

  template <int TDimension>
  static Scalar*
  fghAccessor(CellAccessorS const& accessor) {
    return &accessor.currentCell()->fgh(TDimension);
  }

  static VelocityS*
  velocityAccessor(CellAccessorS const& accessor) {
    return &accessor.currentCell()->velocity();
  }

  static PressureS*
  pressureAccessor(CellAccessorS const& accessor) {
    return &accessor.currentCell()->pressure();
  }

  FluidSimulation::Configuration* _configuration;
  Simulation*                     _simulation;
  GridS*                          _grid;
  ParallelTopologyS*              _parallelDistribution;
  GhostCellsHandlerS*             _handlers;
  VelocityS*                      _maxVelocity;
};
}
}
#endif
