#ifndef FsiSimulation_FluidSimulation_Solvers_LinearSolver_hpp
#define FsiSimulation_FluidSimulation_Solvers_LinearSolver_hpp

#include "GhostLayer/Handlers.hpp"
#include "Grid.hpp"
#include "ParallelDistribution.hpp"
#include "Parameters.hpp"
#include "Private/petscgenerics.hpp"
#include "StructuredMemory/Pointers.hpp"
#include "functions.hpp"
#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

#include <Uni/Logging/macros>

#include <functional>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TGrid,
          typename TStencilGenerator,
          typename TRhsGenerator,
          typename TResultAcquirer>
class LinearSolver {
public:
  typedef TGrid                               GridType;
  typedef typename GridType::CellAccessorType CellAccessorType;
  typedef typename CellAccessorType::CellType CellType;
  typedef typename CellType::Scalar           Scalar;

  enum {
    Dimensions = CellType::Dimensions
  };

  typedef ParallelDistribution<Dimensions>
    ParallelDistributionType;

  typedef Parameters<Scalar, Dimensions> ParametersType;

  typedef
    GhostLayer::LsStencilGenerator::FunctorStack<Dimensions>
    GhostStencilGenerator;
  typedef
    GhostLayer::PetscExchange::FunctorStack<Dimensions>
    GhostRhsGenerator;
  typedef
    GhostLayer::PetscExchange::FunctorStack<Dimensions>
    GhostRhsAcquierer;

  typedef Eigen::Matrix<Scalar, Dimensions, 1> VectorDsType;
  typedef Eigen::Matrix<int, Dimensions, 1>    VectorDiType;

public:
  LinearSolver() {}

  void
  initialize(GridType const*                 grid,
             ParallelDistributionType const* parallelDistribution,
             ParametersType const*           parameters,
             GhostStencilGenerator const*    ghostStencilGenerator,
             GhostRhsGenerator const*        ghostRhsGenerator,
             GhostRhsAcquierer const*        ghostRhsAcquierer,
             Scalar const*                   dt) {
    _grid                  = grid;
    _parallelDistribution  = parallelDistribution;
    _parameters            = parameters;
    _ghostStencilGenerator = ghostStencilGenerator;
    _ghostRhsGenerator     = ghostRhsGenerator;
    _ghostRhsAcquierer     = ghostRhsAcquierer;
    _dt                    = dt;
    KSPCreate(PETSC_COMM_WORLD, &_context);
    PCCreate(PETSC_COMM_WORLD, &_preconditioner);

    VectorDConstPetscIntPointer<Dimensions> localSizes;

    for (int d = 0; d < Dimensions; ++d) {
      auto array = new PetscInt[_parallelDistribution->processorSize(d)];
      localSizes(d) = UniqueConstPetscIntArray(array);

      for (int j = 0; j < _parallelDistribution->processorSize(d); ++j) {
        array[j] = _parallelDistribution->localCellSize(d);
      }
      ++array[0];
      ++array[_parallelDistribution->processorSize(d) - 1];
    }

    DMCreate<Dimensions>(PETSC_COMM_WORLD,
                         createDMBoundaries<Dimensions>(),
                         DMDA_STENCIL_STAR,
                         _parallelDistribution->globalCellSize + 2 *
                         VectorDiType::Ones(),
                         _parallelDistribution->processorSize,
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
    typedef StructuredMemory::Pointers<PetscScalar, Dimensions> Pointers;

    KSPSetComputeRHS(_context, computeRHS, this);
    KSPSolve(_context, PETSC_NULL, _x);

    typename Pointers::Type array;
    DMDAVecGetArray(_da, _x, &array);

    auto corner = _parallelDistribution->corner;

    for (auto const& accessor : _grid->innerGrid) {
      auto index = accessor.indexValues();
      index += corner;

      TResultAcquirer::set(accessor, Pointers::dereference(array, index));
    }

    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        (*_ghostRhsAcquierer)[d][d2](array);
      }
    }

    DMDAVecRestoreArray(_da, _x, &array);
  }

private:
  static PetscErrorCode
  computeMatrix(KSP ksp, Mat A, Mat pc, void* ctx) {
    auto solver = static_cast<LinearSolver*>(ctx);

    PetscScalar stencil[2 * Dimensions + 1];
    MatStencil  row;
    MatStencil  columns[2 * Dimensions + 1];

    for (auto const& accessor : solver->_grid->innerGrid) {
      TStencilGenerator::get(accessor,
                             solver->_parallelDistribution,
                             solver->_parameters,
                             *solver->_dt,
                             stencil,
                             row,
                             columns);

      MatSetValuesStencil(A, 1, &row, 2 * Dimensions + 1, columns, stencil,
                          INSERT_VALUES);
    }

    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        (*solver->_ghostStencilGenerator)[d][d2](A);
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
    typedef StructuredMemory::Pointers<PetscScalar, Dimensions> Pointers;
    typename Pointers::Type array;

    DM da;
    KSPGetDM(ksp, &da);
    DMDAVecGetArray(da, b, &array);

    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        (*solver->_ghostRhsGenerator)[d][d2](array);
      }
    }

    auto corner = solver->_parallelDistribution->corner;

    for (auto const& accessor : solver->_grid->innerGrid) {
      auto value = TRhsGenerator::get(accessor, *solver->_dt);

      auto index = accessor.indexValues();
      index                              += corner;
      Pointers::dereference(array, index) = value;
    }

    DMDAVecRestoreArray(da, b, &array);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    return 0;
  }

  GridType const*                 _grid;
  ParallelDistributionType const* _parallelDistribution;
  ParametersType const*           _parameters;
  GhostStencilGenerator const*    _ghostStencilGenerator;
  GhostRhsGenerator const*        _ghostRhsGenerator;
  GhostRhsAcquierer const*        _ghostRhsAcquierer;
  Scalar const*                   _dt;

  Vec _x;
  DM  _da;
  KSP _context;
  PC  _preconditioner;
};
}
}

#endif
