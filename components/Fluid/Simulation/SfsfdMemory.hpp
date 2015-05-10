#ifndef Fluid_Simulation_SfsfdMemory_hpp
#define Fluid_Simulation_SfsfdMemory_hpp

#include "FsfdMemory.hpp"
#include "SolverTraits.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class SfsfdMemory : public FsfdMemory<TSolverTraits> {
  friend class FsfdMemory<TSolverTraits>;

public:
  using BaseType = FsfdMemory<TSolverTraits>;

  using ScalarType = typename TSolverTraits::ScalarType;

public:
  SfsfdMemory() {}

protected:
  ScalarType _attribute(int const &index,
                        int const &attribute_index,
                        int const &dimension) const {
    return this->BaseType::_attribute(index, attribute_index, dimension);
  }
};
}
}
#endif
