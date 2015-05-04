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
class Grid : public TSolverTraits::SubgridType {
public:
  using SubgridType = typename TSolverTraits::SubgridType;

  using CellAccessorType =  typename TSolverTraits::CellAccessorType;

  using BaseType = SubgridType;

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
  SubgridType                                        innerGrid;
  FactoryType                                        _factory;
  std::array<std::array<SubgridType, 2>, Dimensions> boundaries;
  std::array<std::array<SubgridType, 2>, Dimensions> indentedBoundaries;
};

Fluid_DeclareExternTemplates(Grid);
}
}
#endif
