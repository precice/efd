#ifndef FsiSimulation_Solvers_PoissonSolver_hpp
#define FsiSimulation_Solvers_PoissonSolver_hpp

#include "GhostCellsHandler.hpp"
#include "Grid.hpp"
#include "ParallelTopology.hpp"
#include "PetscStructures.hpp"
#include "StructuredMemory/Pointers.hpp"
#include "stencils/mystencils.hpp"
#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace Solvers {
template <typename TCellAccessor, typename Scalar, int D>
class PoissonSolver {
public:
  typedef Grid<TCellAccessor, D>      SpecializedGrid;
  typedef ParallelTopology<D>         SpecializedParallelTopology;
  typedef GhostCellsHandler<D>        SpecializedGhostCellsHandler;
  typedef Scalar const*               ScalarPointer;
  typedef Eigen::Matrix<Scalar, D, 1> VectorDs;
  typedef Eigen::Matrix<int, D, 1>    VectorDi;

public:
  PoissonSolver() {}

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
        array[j] = _parallelTopology->localSize(d);
      }
      ++array[0];
      ++array[_parallelTopology->processorSize(d) - 1];
    }

    DMCreate<D>(PETSC_COMM_WORLD,
                createDMBoundaries<D>(),
                DMDA_STENCIL_STAR,
                _parallelTopology->globalSize + 2 * VectorDi::Ones(),
                _parallelTopology->processorSize,
                1,
                2,
                localSizes,
                &_da);

    DMCreateGlobalVector(_da, &_x);
    KSPSetDM(_context, _da);
    KSPSetComputeOperators(_context, computeMatrix, this);
    setupCustomOptions(_context, _preconditioner);
  }

  void
  solve() {
    typedef StructuredMemory::Pointers<PetscScalar, D> Pointers;

    KSPSetComputeRHS(_context, computeRHS, this);
    KSPSolve(_context, PETSC_NULL, _x);

    typename Pointers::Type array;
    DMDAVecGetArray(_da, _x, &array);

    for (auto const& accessor : _grid->innerGrid) {
      accessor.currentCell()->pressure() =
        Pointers::dereference(array, accessor.indexValues());
    }

    DMDAVecRestoreArray(_da, _x, &array);
  }

private:
  static PetscErrorCode
  computeMatrix(KSP ksp, Mat A, Mat pc, void* ctx) {
    auto solver = static_cast<PoissonSolver*>(ctx);

    PetscScalar stencil[2 * D + 1];
    MatStencil  row;
    MatStencil  columns[2 * D + 1];

    for (auto const& accessor : solver->_grid->innerGrid) {
      //logInfo("{1}", accessor.indexValues().transpose());
      typedef
        PressurePoissonStencilProcessing<typename SpecializedGrid::Base,
                                         TCellAccessor,
                                         Scalar,
                                         D>
        StencilProcessing;
      StencilProcessing::compute(solver->_grid->innerGrid,
                                 accessor,
                                 stencil,
                                 row,
                                 columns);

      //logInfo("Columns1 {1} {2}", columns[0].i, columns[0].j);
      //logInfo("Columns2 {1} {2}", columns[1].i, columns[1].j);
      //logInfo("Columns3 {1} {2}", columns[2].i, columns[2].j);
      //logInfo("Columns4 {1} {2}", columns[3].i, columns[3].j);
      //logInfo("Columns5 {1} {2}", columns[4].i, columns[4].j);
      //logInfo("Row {1} {2}",      row.i,        row.j);
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
    auto solver = static_cast<PoissonSolver*>(ctx);
    typedef StructuredMemory::Pointers<PetscScalar, D> Pointers;
    typename Pointers::Type array;

    DM da;
    KSPGetDM(ksp, &da);
    DMDAVecGetArray(da, b, &array);

    for (auto const& accessor : solver->_grid->innerGrid) {
      typedef RhsProcessing<TCellAccessor, Scalar, D> Rhs;

      Pointers::dereference(array, accessor.indexValues()) =
        Rhs::compute(accessor, *solver->_dt);
    }

    for (int d = 0; d < D; ++d) {
      //
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

  Vec _x; // ! Petsc vectors for solution and RHS
  DM  _da;                        // ! Topology manager
  KSP _context; // ! Solver context
  PC  _preconditioner;                        // ! Preconditioner

  // Indices for filling the matrices and right hand side
  int _limitsX[2], _limitsY[2], _limitsZ[2];

  PetscInt _firstX, _lengthX, _firstY, _lengthY, _firstZ, _lengthZ;

  // Additional variables used to determine where to write back the results
  int _offsetX, _offsetY, _offsetZ;
};
}
}

#endif
