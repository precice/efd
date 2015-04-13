#ifndef Fluid_Simulation_Grid_hpp
#define Fluid_Simulation_Grid_hpp

#include "IfsfdCellAccessor.hpp"
#include "SfsfdCellAccessor.hpp"
#include "SolverTraits.hpp"

#include <Uni/StructuredGrid/Basic/Grid>

#include <array>
#include <string>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class Grid :
  public Uni::StructuredGrid::Basic::Grid
  <typename TSolverTraits::CellAccessorType> {
public:
  using BaseType = Uni::StructuredGrid::Basic::Grid
                   <typename TSolverTraits::CellAccessorType>;

  using CellAccessorType =  typename TSolverTraits::CellAccessorType;

  enum {
    Dimensions = TSolverTraits::Dimensions
  };

  using VectorDiType = typename BaseType::VectorDi;

  using FactoryType = typename BaseType::Factory;

  using IteratorType = typename BaseType::Iterator;

public:
  Grid() {}

  Grid(Grid const& other) = delete;

  ~Grid() {}

  Grid&
  operator=(Grid const& other) = delete;

  void
  initialize(VectorDiType const&, FactoryType const&);

  std::string
  toString() const;

public:
  BaseType                                        innerGrid;
  FactoryType                                     _factory;
  std::array<std::array<BaseType, 2>, Dimensions> boundaries;
  std::array<std::array<BaseType, 2>, Dimensions> indentedBoundaries;
};

extern template class Grid
  < SfsfdSolverTraits < 0, 0, double, 2 >>;
extern template class Grid
  < SfsfdSolverTraits < 0, 1, double, 2 >>;
extern template class Grid
  < SfsfdSolverTraits < 1, 0, double, 2 >>;
extern template class Grid
  < SfsfdSolverTraits < 1, 1, double, 2 >>;
extern template class Grid
  < SfsfdSolverTraits < 0, 0, double, 3 >>;
extern template class Grid
  < SfsfdSolverTraits < 0, 1, double, 3 >>;
extern template class Grid
  < SfsfdSolverTraits < 1, 0, double, 3 >>;
extern template class Grid
  < SfsfdSolverTraits < 1, 1, double, 3 >>;

extern template class Grid
  < IfsfdSolverTraits < 0, 0, double, 2 >>;
extern template class Grid
  < IfsfdSolverTraits < 0, 1, double, 2 >>;
extern template class Grid
  < IfsfdSolverTraits < 1, 0, double, 2 >>;
extern template class Grid
  < IfsfdSolverTraits < 1, 1, double, 2 >>;
extern template class Grid
  < IfsfdSolverTraits < 0, 0, double, 3 >>;
extern template class Grid
  < IfsfdSolverTraits < 0, 1, double, 3 >>;
extern template class Grid
  < IfsfdSolverTraits < 1, 0, double, 3 >>;
extern template class Grid
  < IfsfdSolverTraits < 1, 1, double, 3 >>;
}
}
#endif
