#ifndef FsiSimulation_Solvers_GhostPressureStencilHandlers_hpp
#define FsiSimulation_Solvers_GhostPressureStencilHandlers_hpp

#include "ParallelTopology.hpp"

#include <petscdm.h>
#include <petscdmda.h>

#include <array>
#include <functional>

namespace FsiSimulation {
namespace Solvers {
template <typename TGrid,
          typename Scalar,
          int D,
          int Dimension,
          int Direction>
class GhostPressureStencilHandler {
public:
  typedef ParallelTopology<D> SpecializedParallelTopology;

public:
  static void
  dirichletPressureStencil(PetscScalar* stencil) {
    stencil[0] = 1;
    stencil[1] = -1;
  }

  static void
  dirichlet(TGrid const*                       grid,
            SpecializedParallelTopology const* parallelTopology,
            Mat&                               A) {
    genericHandler<dirichletPressureStencil>(grid,
                                             parallelTopology,
                                             A);
  }

  template <void(* TStencil) (PetscScalar*)>
  static void
  genericHandler(TGrid const*                       grid,
                 SpecializedParallelTopology const* parallelTopology,
                 Mat&                               A) {
    PetscScalar stencil[2];
    MatStencil  row;
    MatStencil  columns[2];

    for (auto const& accessor : grid->boundaries[Dimension][Direction]) {
      if (Dimension == 0) {
        if (Direction == 1) {
          columns[0].i = parallelTopology->globalSize(0) - 1;
          columns[1].i = columns[0].i - 1;
        } else {
          columns[0].i = 0;
          columns[1].i = columns[0].i + 1;
        }
        columns[0].j = accessor.indexValue(1);
        columns[1].j = columns[0].j;

        if (D == 3) {
          columns[0].k = accessor.indexValue(2);
          columns[1].k = columns[0].k;
        }
      } else if (Dimension == 1) {
        if (Direction == 1) {
          columns[0].j = parallelTopology->globalSize(1) - 1;
          columns[1].j = columns[0].j - 1;
        } else {
          columns[0].j = 0;
          columns[1].j = columns[0].j + 1;
        }
        columns[0].i = accessor.indexValue(0);
        columns[1].i = columns[0].i;

        if (D == 3) {
          columns[0].k = accessor.indexValue(2);
          columns[1].k = columns[0].k;
        }
      } else if (Dimension == 2) {
        if (Direction == 1) {
          columns[0].k = parallelTopology->globalSize(2) - 1;
          columns[1].k = columns[0].k - 1;
        } else {
          columns[0].k = 0;
          columns[1].k = columns[0].k + 1;
        }
        columns[0].i = accessor.indexValue(0);
        columns[1].i = columns[0].i;

        columns[0].j = accessor.indexValue(1);
        columns[1].j = columns[0].j;
      }

      row = columns[0];

      TStencil(stencil);
      MatSetValuesStencil(A, 1, &row, 2, columns, stencil, INSERT_VALUES);
    }
  }
};
template <typename TGrid, typename Scalar, int D>
class GhostPressureStencilHandlers {
public:
  template <int TDimension, int TDirection>
  using PressureStencilHandler = GhostPressureStencilHandler<
          TGrid, Scalar, D, TDimension, TDirection>;
  typedef ParallelTopology<D> SpecializedParallelTopology;
  typedef std::function<void (TGrid const*,
                              SpecializedParallelTopology const*,
                              Mat&)> Handler;
  template <typename THandler>
  using HandlerStack = std::array<std::array<THandler, 2>, D>;
  typedef HandlerStack<Handler> Stack;
};

template <typename TGrid, typename Scalar, int D>
class DirichletStack {};

template <typename TGrid, typename Scalar>
class DirichletStack<TGrid, Scalar, 2> {
public:
  typedef GhostPressureStencilHandlers<TGrid, Scalar, 2> Source;
  typedef typename Source::Stack                         Stack;
  typedef typename Source::Handler                       Handler;

  static Stack const&
  get() {
    static Stack _instance =
    { Handler(Source::template PressureStencilHandler<0, 0>::dirichlet),
      Handler(Source::template PressureStencilHandler<0, 1>::dirichlet),
      Handler(Source::template PressureStencilHandler<1, 0>::dirichlet),
      Handler(Source::template PressureStencilHandler<1, 1>::dirichlet) };

    return _instance;
  }
};

template <typename TGrid, typename Scalar>
class DirichletStack<TGrid, Scalar, 3> {
public:
  typedef GhostPressureStencilHandlers<TGrid, Scalar, 3> Source;
  typedef typename Source::Stack                         Stack;
  typedef typename Source::Handler                       Handler;

  static Stack const&
  get() {
    static Stack _instance =
    { Handler(Source::template PressureStencilHandler<0, 0>::dirichlet),
      Handler(Source::template PressureStencilHandler<0, 1>::dirichlet),
      Handler(Source::template PressureStencilHandler<1, 0>::dirichlet),
      Handler(Source::template PressureStencilHandler<1, 1>::dirichlet),
      Handler(Source::template PressureStencilHandler<2, 0>::dirichlet),
      Handler(Source::template PressureStencilHandler<2, 1>::dirichlet) };

    return _instance;
  }
};
}
}
#endif
