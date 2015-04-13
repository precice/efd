#ifndef Fluid_Simulation_SfsfdMemory_hpp
#define Fluid_Simulation_SfsfdMemory_hpp

#include "FsfdDebugMemory.hpp"
#include "FsfdMemory.hpp"
#include "SolverTraits.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class SfsfdMemory;

template <typename TSolverTraits>
class SfsfdDebugMemory;

template <typename TSolverTraits>
struct SfsfdMemoryTraits {
  using Type = SfsfdMemory<TSolverTraits>;

  enum {
    AdditionalAttributeSize
      = FsfdMemoryTraits<TSolverTraits>::AdditionalAttributeSize + 0
  };
};

template <typename TSolverTraits>
struct SfsfdDebugMemoryTraits {
  using Type = SfsfdDebugMemory<TSolverTraits>;

  enum {
    AdditionalAttributeSize
      = FsfdDebugMemoryTraits<TSolverTraits>::AdditionalAttributeSize + 0
  };
};

template <typename TSolverTraits>
class SfsfdMemory : public FsfdMemory
                    < TSolverTraits,
                           SfsfdMemoryTraits < TSolverTraits >> {
  friend class FsfdMemory
               < TSolverTraits, SfsfdMemoryTraits < TSolverTraits >>;
public:
  using Base = FsfdMemory
               < TSolverTraits, SfsfdMemoryTraits < TSolverTraits >>;

  using ScalarType = typename TSolverTraits::ScalarType;

public:
  SfsfdMemory() {}
};

template <typename TSolverTraits>
class SfsfdDebugMemory : public FsfdDebugMemory
                         < TSolverTraits,
                                SfsfdDebugMemoryTraits < TSolverTraits >> {
  friend class FsfdMemory
               < TSolverTraits, SfsfdDebugMemoryTraits < TSolverTraits >>;

public:
  using Base = FsfdDebugMemory
               < TSolverTraits, SfsfdDebugMemoryTraits < TSolverTraits >>;

  using ScalarType = typename TSolverTraits::ScalarType;

public:
  SfsfdDebugMemory() {}
};

extern template class SfsfdMemory
  < SfsfdSolverTraits < 0, 0, double, 2 >>;
extern template class SfsfdMemory
  < SfsfdSolverTraits < 0, 1, double, 2 >>;
extern template class SfsfdMemory
  < SfsfdSolverTraits < 1, 0, double, 2 >>;
extern template class SfsfdMemory
  < SfsfdSolverTraits < 1, 1, double, 2 >>;
extern template class SfsfdMemory
  < SfsfdSolverTraits < 0, 0, double, 3 >>;
extern template class SfsfdMemory
  < SfsfdSolverTraits < 0, 1, double, 3 >>;
extern template class SfsfdMemory
  < SfsfdSolverTraits < 1, 0, double, 3 >>;
extern template class SfsfdMemory
  < SfsfdSolverTraits < 1, 1, double, 3 >>;
}
}
#endif
