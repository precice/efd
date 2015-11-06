#include "createinctance.hpp"

#include "Writer.hpp"

#include "Simulation/IterationResultWriter.hpp"
#include "Simulation/IfsfdMemory.hpp"
#include "Simulation/SfsfdMemory.hpp"

namespace Fluid {
namespace Simulation {
namespace XdmfHdf5Output {
template <typename TMemory>
std::unique_ptr<IterationResultWriter>
create_instance(TMemory const*                 memory) {
  return std::unique_ptr<IterationResultWriter>(new Writer<TMemory>(memory));
}

template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <0, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <1, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <0, double, 3>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <1, double, 3>
                ::MemoryType const*);

template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <0, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <1, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <0, double, 3>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <1, double, 3>
                ::MemoryType const*);
}
}
}
