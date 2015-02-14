#ifndef FsiSimulation_FluidSimulation_GhostLayer_PeStencilGenerator_Handler_hpp
#define FsiSimulation_FluidSimulation_GhostLayer_PeStencilGenerator_Handler_hpp

#include "Private/utilities.hpp"

#include "FluidSimulation/ParallelDistribution.hpp"

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
  typedef typename GridType::CellAccessor     CellAccessorType;
  typedef typename CellAccessorType::CellType CellType;

  enum {
    Dimensions = CellType::Dimensions
  };

  typedef ParallelDistribution<Dimensions> SpecializedParallelTopology;

public:
  Handler(TGrid const*                       grid,
          SpecializedParallelTopology const* parallelTopology)
    : _grid(grid),
      _parallelTopology(parallelTopology) {}

  ~Handler() {
    logInfo("GhostPeStencilGeneratorHanlers destroyed");
  }

  void
  dirichletLeft(Mat& A) {
    genericHandler<Handler::dirichletStencilLeft>(A);
  }

  void
  dirichletRight(Mat& A) {
    genericHandler<Handler::dirichletStencilRight>(A);
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
    genericHandler<Handler::neumannStencilRight>(A);
  }

  void
  neumannMiddle(Mat& A) {
    genericHandler<Handler::neumannStencilMiddle>(A);
  }

  template <void(* TStencil) (PetscScalar*)>
  void
  genericHandler(Mat& A) {
    PetscScalar stencil[3];
    MatStencil  row;
    MatStencil  columns[3];

    auto corner = _parallelTopology->corner;

    for (auto const
         & accessor : _grid->indentedBoundaries[TDimension][TDirection]) {
      if (TDimension == 0) {
        if (TDirection == 1) {
          columns[0].i = _parallelTopology->globalCellSize(0)
                         + _grid->innerGrid.leftIndent() (0);
          columns[1].i = columns[0].i - 1;
          columns[2].i = columns[0].i - 2;
        } else {
          columns[0].i = _grid->innerGrid.leftIndent() (0) - 1;
          columns[1].i = columns[0].i + 1;
          columns[2].i = columns[0].i + 2;
        }
        columns[0].j = accessor.indexValue(1) + corner(1);
        columns[1].j = columns[0].j;
        columns[2].j = columns[0].j;

        if (Dimensions == 3) {
          columns[0].k = accessor.indexValue(2) + corner(2);
          columns[1].k = columns[0].k;
          columns[2].k = columns[0].k;
        }
      } else if (TDimension == 1) {
        if (TDirection == 1) {
          columns[0].j = _parallelTopology->globalCellSize(1)
                         + _grid->innerGrid.leftIndent() (1);
          columns[1].j = columns[0].j - 1;
          columns[2].j = columns[0].j - 2;
        } else {
          columns[0].j = _grid->innerGrid.leftIndent() (1) - 1;
          columns[1].j = columns[0].j + 1;
          columns[2].j = columns[0].j + 2;
        }
        columns[0].i = accessor.indexValue(0) + corner(0);
        columns[1].i = columns[0].i;
        columns[2].i = columns[0].i;

        if (Dimensions == 3) {
          columns[0].k = accessor.indexValue(2) + corner(2);
          columns[1].k = columns[0].k;
          columns[2].k = columns[0].k;
        }
      } else if (TDimension == 2) {
        if (TDirection == 1) {
          columns[0].k = _parallelTopology->globalCellSize(2)
                         + _grid->innerGrid.leftIndent() (2);
          columns[1].k = columns[0].k - 1;
          columns[2].k = columns[0].k - 2;
        } else {
          columns[0].k = _grid->innerGrid.leftIndent() (2) - 1;
          columns[1].k = columns[0].k + 1;
          columns[2].k = columns[0].k + 2;
        }
        columns[0].i = accessor.indexValue(0) + corner(0);
        columns[1].i = columns[0].i;
        columns[2].i = columns[0].i;

        columns[0].j = accessor.indexValue(1) + corner(1);
        columns[1].j = columns[0].j;
        columns[2].j = columns[0].j;
      }

      row = columns[0];

      TStencil(stencil);

      MatSetValuesStencil(A, 1, &row, 3, columns, stencil, INSERT_VALUES);
    }
  }

  static Functor
  getNeumannLeft(TGrid const*                       grid,
                    SpecializedParallelTopology const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology));

    return Functor(std::bind(&Handler::neumannLeft, pointer, _1));
  }

  static Functor
  getNeumannRight(TGrid const*                       grid,
                    SpecializedParallelTopology const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology));

    return Functor(std::bind(&Handler::neumannRight, pointer, _1));
  }

  static Functor
  getNeumannMiddle(TGrid const*                       grid,
                    SpecializedParallelTopology const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology));

    return Functor(std::bind(&Handler::neumannMiddle, pointer, _1));
  }

  static Functor
  getDirichletLeft(TGrid const*                       grid,
                      SpecializedParallelTopology const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology));

    return Functor(std::bind(&Handler::dirichletLeft, pointer, _1));
  }

  static Functor
  getDirichletRight(TGrid const*                       grid,
                      SpecializedParallelTopology const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology));

    return Functor(std::bind(&Handler::dirichletRight, pointer, _1));
  }

  static Functor
  getDirichletMiddle(TGrid const*                       grid,
                      SpecializedParallelTopology const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology));

    return Functor(std::bind(&Handler::dirichletMiddle, pointer, _1));
  }

  static void
  neumannStencilLeft(PetscScalar* stencil) {
    stencil[0] = 1.0;
    stencil[1] = -1.0;
    stencil[2] = 0.0;
  }

  static void
  neumannStencilRight(PetscScalar* stencil) {
    stencil[0] = 1.0;
    stencil[1] = -1.0;
    stencil[2] = 0.0;
  }

  static void
  neumannStencilMiddle(PetscScalar* stencil) {
    stencil[0] = 1.0;
    stencil[1] = -1.0;
    stencil[2] = 0.0;
  }

  static void
  dirichletStencilLeft(PetscScalar* stencil) {
    stencil[0] = 0.5;
    stencil[0] = 0.5;
    stencil[2] = 0.0;
  }

  static void
  dirichletStencilRight(PetscScalar* stencil) {
    stencil[0] = 0.5;
    stencil[0] = 0.5;
    stencil[2] = 0.0;
  }

  static void
  dirichletStencilMiddle(PetscScalar* stencil) {
    stencil[0] = 0.5;
    stencil[1] = 0.5;
    stencil[2] = 0.0;
  }

private:
  TGrid const*                       _grid;
  SpecializedParallelTopology const* _parallelTopology;
};
}
}
}
}
#endif
