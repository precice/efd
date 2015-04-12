#ifndef Fluid_Simulation_IfsfdMemory_hpp
#define Fluid_Simulation_IfsfdMemory_hpp

#include "FsfdDebugMemory.hpp"
#include "FsfdMemory.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits, typename TBaseMemoryTraits>
struct BaseIfsfdMemoryTraits {
  enum {
    AdditionalAttributeSize
      = TBaseMemoryTraits::AdditionalAttributeSize + 0
  };
};

template <typename TSolverTraits, typename TBaseMemory>
class BaseIfsfdMemory : public TBaseMemory {
public:
  using BaseType = TBaseMemory;

  enum {
    Dimensions = TSolverTraits::Dimensions,
    AttributeSize = TBaseMemory::AttributeSize + 0
  };

  using GridGeometryType = typename TSolverTraits::GridGeometryType;

  using ParallelDistributionType
          = typename TSolverTraits::ParallelDistributionType;

  using ParametersType = typename TSolverTraits::ParametersType;

  using GridType = typename TSolverTraits::GridType;

  using BaseGridType = typename TSolverTraits::BaseGridType;

  using CellAccessorType = typename TSolverTraits::CellAccessorType;

  using VectorDsType = typename TSolverTraits::VectorDsType;

  using VectorDiType = typename TSolverTraits::VectorDiType;

  using ScalarType = typename TSolverTraits::ScalarType;

public:
  BaseIfsfdMemory() {}

  void
  initialize(VectorDiType const& processor_size,
             VectorDiType const& global_cell_size,
             VectorDsType const& geometry_width) {
    this->BaseType::initialize(processor_size,
                               global_cell_size,
                               geometry_width);

    _projectionTerm.reset(new ScalarType[this->grid()->size().prod()]);
  }

  ScalarType const*
  projectionTerm() const {
    return _projectionTerm.get();
  }

  ScalarType*
  projectionTerm() {
    return _projectionTerm.get();
  }

  ScalarType const&
  projectionTerm(std::size_t const& index) const {
    return _projectionTerm.get()[index];
  }

  ScalarType&
  projectionTerm(std::size_t const& index) {
    return _projectionTerm.get()[index];
  }

private:
  std::unique_ptr<ScalarType[]> _projectionTerm;
};

template <typename TSolverTraits>
class IfsfdMemory;

template <typename TSolverTraits>
class IfsfdDebugMemory;

template <typename TSolverTraits>
struct IfsfdMemoryTraits {
  using Type = IfsfdMemory<TSolverTraits>;

  enum {
    AdditionalAttributeSize
      = BaseIfsfdMemoryTraits
        < TSolverTraits,
    FsfdMemoryTraits < TSolverTraits >> ::AdditionalAttributeSize + 0
  };
};

template <typename TSolverTraits>
struct IfsfdDebugMemoryTraits {
  using Type = IfsfdDebugMemory<TSolverTraits>;

  enum {
    AdditionalAttributeSize
      = BaseIfsfdMemoryTraits
        < TSolverTraits,
    FsfdDebugMemoryTraits < TSolverTraits >> ::AdditionalAttributeSize + 0
  };
};

template <typename TSolverTraits>
class IfsfdMemory :
  public BaseIfsfdMemory
  < TSolverTraits,
         FsfdMemory < TSolverTraits,
         IfsfdMemoryTraits<TSolverTraits >>> {
  friend class FsfdMemory < TSolverTraits,
    IfsfdMemoryTraits < TSolverTraits >>;
};

template <typename TSolverTraits>
class IfsfdDebugMemory :
  public BaseIfsfdMemory
  < TSolverTraits,
         FsfdDebugMemory < TSolverTraits,
         IfsfdDebugMemoryTraits<TSolverTraits >>> {
  friend class FsfdMemory < TSolverTraits,
    IfsfdDebugMemoryTraits < TSolverTraits >>;
};
}
}
#endif
