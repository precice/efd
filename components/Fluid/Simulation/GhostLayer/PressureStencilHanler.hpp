#pragma once

#include "Private/utilities.hpp"

#include "Simulation/ParallelDistribution.hpp"

#include <Uni/Logging/macros>

#include <petscdm.h>
#include <petscdmda.h>

#include <functional>
#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
namespace LsStencilGenerator {
typedef std::function<void (Mat&)> Functor;
template <int Dimensions>
using FunctorStack = FunctorStack<Functor, Dimensions>;

template <int Dimensions>
Functor
getEmptyFunctor() {
  return Functor([] (Mat&) {});
}

template <typename TGrid,
          int TDimension,
          int TDirection>
class Handler {
public:
  typedef TGrid                               GridType;
  typedef typename GridType::CellAccessorType CellAccessorType;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

  typedef ParallelDistribution<Dimensions> ParallelDistributionType;

public:
  Handler(TGrid const*                    grid,
          ParallelDistributionType const* parallelTopology,
          int const&                      offset = 0)
    : _grid(grid),
    _parallelDistribution(parallelTopology),
    _offset(offset) {}

  ~Handler() {}

  void
  dirichletLeft(Mat& A) {
    genericHandler<Handler::dirichletStencilLeft>(A);
  }

  void
  dirichletRight(Mat& A) {
    genericHandler<Handler::dirichletStencilLeft>(A);
  }

  void
  dirichletMiddle(Mat& A) {
    genericHandler<Handler::dirichletStencilMiddle>(A);
  }

  void
  neumannLeft(Mat& A) {
    genericHandler<Handler::neumannStencilLeft>(A);
  }

  void
  neumannRight(Mat& A) {
    genericHandler<Handler::neumannStencilLeft>(A);
  }

  void
  neumannMiddle(Mat& A) {
    genericHandler<Handler::neumannStencilMiddle>(A);
  }

  template <typename TCellAccessor>
  inline void
  computeGhostIndexes(TCellAccessor const& accessor,
                      MatStencil*          columns,
                      MatStencil*          row) {
    auto corner = _parallelDistribution->corner;

    if (TDimension == 0) {
      if (TDirection == 1) {
        columns[0].i = _parallelDistribution->globalCellSize(0)
                       + _grid->innerGrid.leftIndent() (0);
        columns[1].i = columns[0].i - 1;
      } else {
        columns[0].i = _grid->innerGrid.leftIndent() (0) - 1;
        columns[1].i = columns[0].i + 1;
      }
      columns[0].j = accessor.indexValue(1) + corner(1);
      columns[1].j = columns[0].j;

      if (Dimensions == 3) {
        columns[0].k = accessor.indexValue(2) + corner(2);
        columns[1].k = columns[0].k;
      }
    } else if (TDimension == 1) {
      if (TDirection == 1) {
        columns[0].j = _parallelDistribution->globalCellSize(1)
                       + _grid->innerGrid.leftIndent() (1);
        columns[1].j = columns[0].j - 1;
      } else {
        columns[0].j = _grid->innerGrid.leftIndent() (1) - 1;
        columns[1].j = columns[0].j + 1;
      }
      columns[0].i = accessor.indexValue(0) + corner(0);
      columns[1].i = columns[0].i;

      if (Dimensions == 3) {
        columns[0].k = accessor.indexValue(2) + corner(2);
        columns[1].k = columns[0].k;
      }
    } else if (TDimension == 2) {
      if (TDirection == 1) {
        columns[0].k = _parallelDistribution->globalCellSize(2)
                       + _grid->innerGrid.leftIndent() (2);
        columns[1].k = columns[0].k - 1;
      } else {
        columns[0].k = _grid->innerGrid.leftIndent() (2) - 1;
        columns[1].k = columns[0].k + 1;
      }
      columns[0].i = accessor.indexValue(0) + corner(0);
      columns[1].i = columns[0].i;

      columns[0].j = accessor.indexValue(1) + corner(1);
      columns[1].j = columns[0].j;
    }

    *row = columns[0];
  }

  template <typename TCellAccessor>
  inline void
  computeInnerIndexes(TCellAccessor const& accessor,
                      MatStencil*          columns,
                      MatStencil*          row,
                      PetscScalar*         stencil) {
    auto corner = _parallelDistribution->corner;

    int offset = 0;

    for (int d = 0; d < Dimensions; ++d) {
      typename CellAccessorType::VectorDiType leftIndex;
      typename CellAccessorType::VectorDiType rightIndex;
      int currentOffset = 1;

      if (TDimension == d) {
        if (TDirection == 1) {
          leftIndex  = (accessor.relativeIndex(d, -1) + corner).eval();
          rightIndex = (accessor.relativeIndex(d, +1) + corner).eval();
        } else {
          leftIndex  = (accessor.relativeIndex(d, +1) + corner).eval();
          rightIndex = (accessor.relativeIndex(d, -1) + corner).eval();
        }
      } else {
        currentOffset += 2 + 2 * offset;
        ++offset;
        leftIndex  = (accessor.relativeIndex(d, -1) + corner).eval();
        rightIndex = (accessor.relativeIndex(d, +1) + corner).eval();
      }

      if (TDirection == 1) {
        leftIndex(TDimension)  -= _offset;
        rightIndex(TDimension) -= _offset;
      }

      columns[currentOffset].i     = leftIndex(0);
      columns[currentOffset].j     = leftIndex(1);
      columns[currentOffset + 1].i = rightIndex(0);
      columns[currentOffset + 1].j = rightIndex(1);

      if (Dimensions == 3) {
        columns[currentOffset].k     = leftIndex(2);
        columns[currentOffset + 1].k = rightIndex(2);
      }
      stencil[currentOffset]     = 0.0;
      stencil[currentOffset + 1] = 0.0;
    }
    auto currentIndex = (accessor.index() + corner).eval();

    if (TDirection == 1) {
      currentIndex(TDimension) -= _offset;
    }
    columns[0].i = currentIndex(0);
    columns[0].j = currentIndex(1);

    if (Dimensions == 3) {
      columns[0].k = currentIndex(2);
    }
    *row = columns[0];
  }

  template <void(* TStencil) (PetscScalar*)>
  void
  genericHandler(Mat& A) {
    if (_offset == 0) {
      PetscScalar stencil[2];
      MatStencil  columns[2];
      MatStencil  row;

      TStencil(stencil);

      for (auto const& accessor :
           _grid->indentedBoundaries[TDimension][TDirection]) {
        computeGhostIndexes(accessor, columns, &row);

        MatSetValuesStencil(A, 1, &row, 2, columns, stencil, INSERT_VALUES);
      }
    } else {
      auto corner = _parallelDistribution->corner;

      PetscScalar const_stencil = 1.0;
      MatStencil  const_column;

      PetscScalar stencil[2 * Dimensions + 1];
      MatStencil  columns[2 * Dimensions + 1];
      MatStencil  row;

      for (auto const
           & accessor : _grid->indentedBoundaries[TDimension][TDirection]) {
        computeInnerIndexes(accessor, columns, &row, stencil);
        TStencil(stencil);
        MatSetValuesStencil(A, 1,
                            &row,
                            2 * Dimensions + 1,
                            columns,
                            stencil,
                            INSERT_VALUES);

        const_column.i = accessor.indexValue(0) + corner(0);
        const_column.j = accessor.indexValue(1) + corner(1);

        if (Dimensions == 3) {
          const_column.k = accessor.indexValue(2) + corner(2);
        }
        MatSetValuesStencil(A, 1,
                            &const_column,
                            1,
                            &const_column,
                            &const_stencil,
                            INSERT_VALUES);
      }
    }
  }

  static Functor
  getNeumannLeft(TGrid const*                    grid,
                 ParallelDistributionType const* parallelTopology,
                 int const&                      offset = 0) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology, offset));

    return Functor(std::bind(&Handler::neumannLeft, pointer, _1));
  }

  static Functor
  getNeumannRight(TGrid const*                    grid,
                  ParallelDistributionType const* parallelTopology,
                  int const&                      offset = 0) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology, offset));

    return Functor(std::bind(&Handler::neumannLeft, pointer, _1));
  }

  static Functor
  getNeumannMiddle(TGrid const*                    grid,
                   ParallelDistributionType const* parallelTopology,
                   int const&                      offset = 0) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology, offset));

    return Functor(std::bind(&Handler::neumannMiddle, pointer, _1));
  }

  static Functor
  getDirichletLeft(TGrid const*                    grid,
                   ParallelDistributionType const* parallelTopology,
                   int const&                      offset = 0) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology, offset));

    return Functor(std::bind(&Handler::dirichletLeft, pointer, _1));
  }

  static Functor
  getDirichletRight(TGrid const*                    grid,
                    ParallelDistributionType const* parallelTopology,
                    int const&                      offset = 0) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology, offset));

    return Functor(std::bind(&Handler::dirichletLeft, pointer, _1));
  }

  static Functor
  getDirichletMiddle(TGrid const*                    grid,
                     ParallelDistributionType const* parallelTopology,
                     int const&                      offset = 0) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology, offset));

    return Functor(std::bind(&Handler::dirichletMiddle, pointer, _1));
  }

  static void
  neumannStencilLeft(PetscScalar* stencil) {
    stencil[0] = 1.0;
    stencil[1] = -1.0;
  }

  static void
  neumannStencilMiddle(PetscScalar* stencil) {
    stencil[0] = 1.0;
    stencil[1] = -1.0;
  }

  static void
  dirichletStencilLeft(PetscScalar* stencil) {
    stencil[0] = 1.0;
    stencil[1] = 0.0;
  }

  static void
  dirichletStencilMiddle(PetscScalar* stencil) {
    stencil[0] = 0.5;
    stencil[1] = 0.5;
  }

private:
  TGrid const*                    _grid;
  ParallelDistributionType const* _parallelDistribution;
  int const                       _offset;
};
}
}
}
}
