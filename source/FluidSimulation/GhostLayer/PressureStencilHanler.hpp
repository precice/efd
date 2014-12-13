#ifndef FsiSimulation_Solvers_GhostPressureStencil_Handlers_hpp
#define FsiSimulation_Solvers_GhostPressureStencil_Handlers_hpp

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
namespace PressureStencil {
typedef std::function<void (Mat&)> Functor;
template <int D>
using FunctorStack = FunctorStack<Functor, D>;

template <int D>
Functor
getEmptyFunctor() {
  return Functor([] (Mat&) {});
}

template <typename TGrid,
          typename Scalar,
          int D,
          int TDimension,
          int TDirection>
class Handler {
public:
  typedef ParallelDistribution<D> SpecializedParallelTopology;

public:
  Handler(TGrid const*                       grid,
          SpecializedParallelTopology const* parallelTopology)
    : _grid(grid),
      _parallelTopology(parallelTopology) {}

  ~Handler() {
    logInfo("GhostPressureStencilHanlers destroyed");
  }

  static Functor
  getDirichletHandler(TGrid const*                       grid,
                      SpecializedParallelTopology const* parallelTopology) {
    using std::placeholders::_1;

    auto pointer = std::shared_ptr<Handler>
                     (new Handler(grid, parallelTopology));

    return Functor(std::bind(&Handler::dirichlet,
                             pointer,
                             _1));
  }

  void
  dirichlet(Mat& A) {
    genericHandler<Handler::dirichletPressureStencil>(A);
  }

  template <void(* TStencil) (PetscScalar*)>
  void
  genericHandler(Mat& A) {
    PetscScalar stencil[2];
    MatStencil  row;
    MatStencil  columns[2];

    auto corner = _parallelTopology->corner;

    for (auto const& accessor : _grid->indentedBoundaries[TDimension][TDirection]) {
      if (TDimension == 0) {
        if (TDirection == 1) {
          columns[0].i = _parallelTopology->globalCellSize(0) + 1;
          columns[1].i = columns[0].i - 1;
        } else {
          columns[0].i = 0;
          columns[1].i = columns[0].i + 1;
        }
        columns[0].j = accessor.indexValue(1) + corner(1);
        columns[1].j = columns[0].j;

        if (D == 3) {
          columns[0].k = accessor.indexValue(2) + corner(2);
          columns[1].k = columns[0].k;
        }
      } else if (TDimension == 1) {
        if (TDirection == 1) {
          columns[0].j = _parallelTopology->globalCellSize(1) + 1;
          columns[1].j = columns[0].j - 1;
        } else {
          columns[0].j = 0;
          columns[1].j = columns[0].j + 1;
        }
        columns[0].i = accessor.indexValue(0) + corner(0);
        columns[1].i = columns[0].i;

        if (D == 3) {
          columns[0].k = accessor.indexValue(2) + corner(2);
          columns[1].k = columns[0].k;
        }
      } else if (TDimension == 2) {
        if (TDirection == 1) {
          columns[0].k = _parallelTopology->globalCellSize(2) + 1;
          columns[1].k = columns[0].k - 1;
        } else {
          columns[0].k = 0;
          columns[1].k = columns[0].k + 1;
        }
        columns[0].i = accessor.indexValue(0) + corner(0);
        columns[1].i = columns[0].i;

        columns[0].j = accessor.indexValue(1) + corner(1);
        columns[1].j = columns[0].j;
      }

      row = columns[0];

      TStencil(stencil);

      MatSetValuesStencil(A, 1, &row, 2, columns, stencil, INSERT_VALUES);
    }
  }

  static void
  dirichletPressureStencil(PetscScalar* stencil) {
    stencil[0] = 1;
    stencil[1] = -1;
  }

private:
  TGrid const*                       _grid;
  SpecializedParallelTopology const* _parallelTopology;
};

template <typename TGrid, typename Scalar, int D>
class DirichletStack {};

template <typename TGrid, typename Scalar>
class DirichletStack<TGrid, Scalar, 2> {
public:
  typedef ParallelDistribution<2> SpecializedParallelTopology;
  template <int TDimension, int TDirection>
  using _Handler = Handler<TGrid, Scalar, 2, TDimension, TDirection>;

  static FunctorStack<2>
  create(TGrid const*                       grid,
         SpecializedParallelTopology const* topology) {
    FunctorStack<2> _instance =
    { _Handler<0, 0>::getDirichletHandler(grid, topology),
      _Handler<0, 1>::getDirichletHandler(grid, topology),
      _Handler<1, 0>::getDirichletHandler(grid, topology),
      _Handler<1, 1>::getDirichletHandler(grid, topology) };

    return _instance;
  }
};

template <typename TGrid, typename Scalar>
class DirichletStack<TGrid, Scalar, 3> {
public:
  typedef ParallelDistribution<3> SpecializedParallelTopology;
  template <int TDimension, int TDirection>
  using _Handler = Handler<TGrid, Scalar, 3, TDimension, TDirection>;

  static FunctorStack<3>
  create(TGrid const*                       grid,
         SpecializedParallelTopology const* topology) {
    FunctorStack<3> _instance =
    { _Handler<0, 0>::getDirichletHandler(grid, topology),
      _Handler<0, 1>::getDirichletHandler(grid, topology),
      _Handler<1, 0>::getDirichletHandler(grid, topology),
      _Handler<1, 1>::getDirichletHandler(grid, topology),
      _Handler<2, 0>::getDirichletHandler(grid, topology),
      _Handler<2, 1>::getDirichletHandler(grid, topology) };

    return _instance;
  }
};
}
}
}
}
#endif
