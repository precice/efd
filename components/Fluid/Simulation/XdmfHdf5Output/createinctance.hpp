#pragma once

#include "Simulation/SolverTraits.hpp"
#include "Simulation/IterationResultWriter.hpp"

#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
namespace XdmfHdf5Output {
template <typename TMemory>
std::unique_ptr<IterationResultWriter>
create_instance(TMemory const*);

extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 0, 0, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 0, 1, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 1, 0, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 1, 1, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 0, 0, double, 3>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 0, 1, double, 3>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 1, 0, double, 3>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 1, 1, double, 3>
                ::MemoryType const*);

extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 0, 0, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 0, 1, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 1, 0, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 1, 1, double, 2>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 0, 0, double, 3>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 0, 1, double, 3>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 1, 0, double, 3>
                ::MemoryType const*);
extern template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 1, 1, double, 3>
                ::MemoryType const*);
}
}
}
