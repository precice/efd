#pragma once

#include "Simulation/IterationResultWriter.hpp"
#include "Simulation/SolverTraits.hpp"

#include <memory>

namespace Fluid {
namespace Simulation {
namespace VtkOutput {
template <typename TMemory>
std::unique_ptr<IterationResultWriter>
create_instance(TMemory const*);

extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <0, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <1, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <0, double, 3>
                ::MemoryType const*);

extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <0, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <1, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <0, double, 3>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <1, double, 3>
                ::MemoryType const*);
}
}
}
