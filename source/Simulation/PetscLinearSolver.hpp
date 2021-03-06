#ifndef Fluid_Simulation_LinearSolver_hpp
#define Fluid_Simulation_LinearSolver_hpp

#include "LinearSolver.hpp"

#include "ParallelDistribution.hpp"
#include "Private/petscgenerics.hpp"
#include "functions.hpp"

#include "GhostLayer/PetscExchangeHandler.hpp"
#include "GhostLayer/PressureStencilHanler.hpp"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

#include <Uni/Logging/macros>
#include <Uni/StructuredMemory/Pointers>

#include <functional>

namespace Fluid {
namespace Simulation {
template <typename TGrid,
          typename TStencilGenerator,
          typename TRhsGenerator,
          typename TResultAcquirer>
class PetscLinearSolver : public LinearSolver {
public:
  typedef TGrid GridType;

  typedef typename GridType::CellAccessorType CellAccessorType;

  typedef typename CellAccessorType::ScalarType ScalarType;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

  typedef ParallelDistribution<Dimensions> ParallelDistributionType;

  typedef
    GhostLayer::LsStencilGenerator::FunctorStack<Dimensions>
    GhostStencilGenerator;

  typedef
    GhostLayer::PetscExchange::FunctorStack<Dimensions>
    GhostRhsGenerator;

  typedef
    GhostLayer::PetscExchange::FunctorStack<Dimensions>
    GhostRhsAcquierer;

  typedef Eigen::Matrix<ScalarType, Dimensions, 1> VectorDsType;

  typedef Eigen::Matrix<int, Dimensions, 1> VectorDiType;

public:
  PetscLinearSolver() {}

  void
  initialize(GridType const*                 grid,
             ParallelDistributionType const* parallelDistribution,
             GhostStencilGenerator const*    ghostStencilGenerator,
             GhostRhsGenerator const*        ghostRhsGenerator,
             GhostRhsAcquierer const*        ghostRhsAcquierer) {
    _grid                  = grid;
    _parallelDistribution  = parallelDistribution;
    _ghostStencilGenerator = ghostStencilGenerator;
    _ghostRhsGenerator     = ghostRhsGenerator;
    _ghostRhsAcquierer     = ghostRhsAcquierer;
    KSPCreate(PETSC_COMM_WORLD, &_context);
    PCCreate(PETSC_COMM_WORLD, &_preconditioner);

    VectorDConstPetscIntPointer<Dimensions> localSizes;

    for (int d = 0; d < Dimensions; ++d) {
      auto array = new PetscInt[_parallelDistribution->processorSize(d)];
      localSizes(d) = UniqueConstPetscIntArray(array);

      for (int j = 0; j < _parallelDistribution->processorSize(d) - 1; ++j) {
        array[j] = _parallelDistribution->uniformLocalCellSize(d);
      }
      array[_parallelDistribution->processorSize(d) - 1]
        = _parallelDistribution->lastLocalCellSize(d) + 1;
      ++array[0];
    }

    DMCreate<Dimensions>(PETSC_COMM_WORLD,
                         createDMBoundaries<Dimensions>(),
                         DMDA_STENCIL_STAR,
                         _parallelDistribution->globalCellSize
                         + 2 * VectorDiType::Ones(),
                         _parallelDistribution->processorSize,
                         1,
                         2,
                         localSizes,
                         &_da);

    DMCreateGlobalVector(_da, &_x);
    KSPSetDM(_context, _da);
    KSPSetFromOptions(_context);

    KSPGetTolerances(_context,
                     &_rtol,
                     &_atol,
                     &_dtol,
                     &_maxits);
  }

  long double
  absoluteTolerance() const {
    return _atol;
  }

  void
  absoluteTolerance(long double const& tolerance) {
    _atol = tolerance;
  }

  long double
  relativeTolerance() const {
    return _rtol;
  }

  void
  relativeTolerance(long double const& tolerance) {
    _rtol = tolerance;
  }

  void
  update() {
    KSPSetTolerances(_context,
                     _rtol,
                     _atol,
                     _dtol,
                     _maxits);
    KSPSetComputeOperators(_context, computeMatrix, this);
    KSPSetUp(_context);
  }

  void
  solve() {
    using Pointers = Uni::StructuredMemory::Pointers<PetscScalar, Dimensions>;

    KSPSetComputeRHS(_context, computeRHS, this);
    KSPSolve(_context, PETSC_NULL, _x);

    typename Pointers::Type array;
    DMDAVecGetArray(_da, _x, &array);

    auto corner = _parallelDistribution->corner;

    for (auto const& accessor : _grid->innerGrid) {
      auto index = accessor.indexValues();
      index += corner;

      resultAcquierer.set(accessor,
                          Pointers::dereference(array, index));
    }

    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        (*_ghostRhsAcquierer)[d][d2](array);
      }
    }

    DMDAVecRestoreArray(_da, _x, &array);
  }

  void
  release() {
    KSPDestroy(&_context);
  }

  TStencilGenerator stencilGenerator;
  TRhsGenerator     rhsGenerator;
  TResultAcquirer   resultAcquierer;

private:
  static PetscErrorCode
#if ((PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 5))
  computeMatrix(KSP ksp, Mat A, Mat pc, void* ctx) {
#else
  computeMatrix(KSP ksp, Mat A, Mat pc, MatStructure * matStructure,
                void* ctx) {
    ((void)matStructure);
#endif
    ((void)ksp);
    ((void)pc);
    auto this_solver = static_cast<PetscLinearSolver*>(ctx);

    PetscScalar stencil[2 * Dimensions + 1];
    MatStencil  row;
    MatStencil  columns[2 * Dimensions + 1];

    for (auto const& accessor : this_solver->_grid->innerGrid) {
      this_solver->stencilGenerator.get(accessor,
                                        stencil,
                                        row,
                                        columns);

      MatSetValuesStencil(A, 1, &row, 2 * Dimensions + 1, columns, stencil,
                          INSERT_VALUES);
    }

    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        (*this_solver->_ghostStencilGenerator)[d][d2](A);
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
    auto this_solver = static_cast<PetscLinearSolver*>(ctx);
    using Pointers = Uni::StructuredMemory::Pointers<PetscScalar, Dimensions>;
    typename Pointers::Type array;

    DM da;
    KSPGetDM(ksp, &da);
    DMDAVecGetArray(da, b, &array);

    auto corner = this_solver->_parallelDistribution->corner;

    for (auto const& accessor : this_solver->_grid->innerGrid) {
      auto value = this_solver->rhsGenerator.get(accessor);

      auto index = accessor.indexValues();
      index                              += corner;
      Pointers::dereference(array, index) = value;
    }

    for (int d = 0; d < Dimensions; ++d) {
      for (int d2 = 0; d2 < 2; ++d2) {
        (*this_solver->_ghostRhsGenerator)[d][d2](array);
      }
    }

    DMDAVecRestoreArray(da, b, &array);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    return 0;
  }

  GridType const*                 _grid;
  ParallelDistributionType const* _parallelDistribution;
  GhostStencilGenerator const*    _ghostStencilGenerator;
  GhostRhsGenerator const*        _ghostRhsGenerator;
  GhostRhsAcquierer const*        _ghostRhsAcquierer;

  Vec _x;
  DM  _da;
  KSP _context;
  PC  _preconditioner;
  PetscReal _atol;
  PetscReal _rtol;
  PetscReal _dtol;
  PetscInt _maxits;
};
}
}
#endif
