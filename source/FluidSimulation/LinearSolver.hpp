#ifndef FsiSimulation_Solvers_PoissonSolver_hpp
#define FsiSimulation_Solvers_PoissonSolver_hpp

#include "GhostLayer/Handlers.hpp"
#include "Grid.hpp"
#include "ParallelDistribution.hpp"
#include "Private/petscgenerics.hpp"
#include "StructuredMemory/Pointers.hpp"
#include "functions.hpp"
#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TCellAccessor, typename Scalar, int D>
class LinearSolver {
public:
  typedef Grid<TCellAccessor, D>      SpecializedGrid;
  typedef ParallelDistribution<D>         SpecializedParallelTopology;
  typedef typename GhostLayer::Handlers<D>        SpecializedGhostCellsHandler;
  typedef Scalar const*               ScalarPointer;
  typedef Eigen::Matrix<Scalar, D, 1> VectorDs;
  typedef Eigen::Matrix<int, D, 1>    VectorDi;

public:
  LinearSolver() {}

  void
  initialize(SpecializedGrid const*              grid,
             SpecializedParallelTopology const*  parallelTopology,
             SpecializedGhostCellsHandler const* ghostCellsHandler,
             Scalar const*                       dt) {
    _grid              = grid;
    _parallelTopology  = parallelTopology;
    _ghostCellsHandler = ghostCellsHandler;
    _dt                = dt;
    KSPCreate(PETSC_COMM_WORLD, &_context);
    PCCreate(PETSC_COMM_WORLD, &_preconditioner);

    VectorDConstPetscIntPointer<D> localSizes;

    for (int d = 0; d < D; ++d) {
      auto array = new PetscInt[_parallelTopology->processorSize(d)];
      localSizes(d) = UniqueConstPetscIntArray(array);

      for (int j = 0; j < _parallelTopology->processorSize(d); ++j) {
        array[j] = _parallelTopology->localCellSize(d);
      }
      ++array[0];
      ++array[_parallelTopology->processorSize(d) - 1];
    }

    DMCreate<D>(PETSC_COMM_WORLD,
                createDMBoundaries<D>(),
                DMDA_STENCIL_STAR,
                _parallelTopology->globalCellSize + 2 * VectorDi::Ones(),
                _parallelTopology->processorSize,
                1,
                2,
                localSizes,
                &_da);

    DMCreateGlobalVector(_da, &_x);
    KSPSetDM(_context, _da);
    KSPSetComputeOperators(_context, computeMatrix, this);
    setupCustomOptions(_context, _preconditioner);
    DMDALocalInfo petscInfo;
    DMDAGetLocalInfo(_da, &petscInfo);
    logInfo("dim {1} dof {2} sw {3} \n"
            "global number of grid points in each direction {4} {5} {6} \n"
            "starting point of this processor, excluding ghosts {7} {8} {9} \n"
            "number of grid points on this processor, excluding ghosts {10} {11} {12} \n"
            "starting point of this processor including ghosts {13} {14} {15} \n"
            "number of grid points on this processor including ghosts {16} {17} {18} \n",
            petscInfo.dim, petscInfo.dof, petscInfo.sw,
            petscInfo.mx, petscInfo.my, petscInfo.mz,
            petscInfo.xs, petscInfo.ys, petscInfo.zs,
            petscInfo.xm, petscInfo.ym, petscInfo.zm,
            petscInfo.gxs, petscInfo.gys, petscInfo.gzs,
            petscInfo.gxm, petscInfo.gym, petscInfo.gzm);

    int _firstX;
    int _firstY;
    int _firstZ;
    int _lengthX;
    int _lengthY;
    int _lengthZ;
    DMDAGetCorners(_da,
                   &_firstX,
                   &_firstY,
                   &_firstZ,
                   &_lengthX,
                   &_lengthY,
                   &_lengthZ);

    logInfo("Corner beginning {1} {2} {3}\n"
            "Corener end      {4} {5} {6}\n",
            _firstX, _firstY, _firstZ,
            _firstX + _lengthX, _firstY + _lengthY, _firstZ + _lengthZ);
  }

  void
  solve() {
    typedef StructuredMemory::Pointers<PetscScalar, D> Pointers;

    KSPSetComputeRHS(_context, computeRHS, this);
    KSPSolve(_context, PETSC_NULL, _x);

    typename Pointers::Type array;
    DMDAVecGetArray(_da, _x, &array);

    auto corner = _parallelTopology->corner;

    for (auto const& accessor : _grid->innerGrid) {
      auto index = accessor.indexValues();
      index += corner;

      accessor.currentCell()->pressure() = Pointers::dereference(array, index);
    }

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        _ghostCellsHandler->_pressureInitialization[d][d2](array);
      }
    }

    DMDAVecRestoreArray(_da, _x, &array);
  }

private:
  static PetscErrorCode
  computeMatrix(KSP ksp, Mat A, Mat pc, void* ctx) {
    auto solver = static_cast<LinearSolver*>(ctx);

    PetscScalar stencil[2 * D + 1];
    MatStencil  row;
    MatStencil  columns[2 * D + 1];

    for (auto const& accessor : solver->_grid->innerGrid) {
      typedef
        PressurePoissonStencilProcessing<SpecializedParallelTopology,
                                         TCellAccessor,
                                         Scalar,
                                         D>
        StencilProcessing;

      StencilProcessing::compute(solver->_parallelTopology,
                                 accessor,
                                 stencil,
                                 row,
                                 columns);

      MatSetValuesStencil(A, 1, &row, 2 * D + 1, columns, stencil,
                          INSERT_VALUES);
    }

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        solver->_ghostCellsHandler->_pressureStencilStack[d][d2](A);
      }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    MatNullSpace nullspace;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
    MatSetNullSpace(A, nullspace);
    MatNullSpaceDestroy(&nullspace);

    return 0;
  }

  static PetscErrorCode
  computeRHS(KSP ksp, Vec b, void* ctx) {
    auto solver = static_cast<LinearSolver*>(ctx);
    typedef StructuredMemory::Pointers<PetscScalar, D> Pointers;
    typename Pointers::Type array;

    DM da;
    KSPGetDM(ksp, &da);
    DMDAVecGetArray(da, b, &array);

    for (int d = 0; d < D; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        solver->_ghostCellsHandler->_rhsInitialization[d][d2](array);
      }
    }

    auto corner = solver->_parallelTopology->corner;

    for (auto const& accessor : solver->_grid->innerGrid) {
      typedef RhsProcessing<TCellAccessor, Scalar, D> Rhs;

      auto index = accessor.indexValues();
      index += corner;

      auto value = Rhs::compute(accessor, *solver->_dt);
      Pointers::dereference(array, index) = value;
    }

    DMDAVecRestoreArray(da, b, &array);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    return 0;
  }

  SpecializedGrid const*              _grid;
  SpecializedParallelTopology const*  _parallelTopology;
  SpecializedGhostCellsHandler const* _ghostCellsHandler;
  ScalarPointer                       _dt;

  Vec _x;
  DM  _da;
  KSP _context;
  PC  _preconditioner;
};
}
}

#endif
