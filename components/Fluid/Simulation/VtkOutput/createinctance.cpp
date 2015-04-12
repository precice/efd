#include "createinctance.hpp"

#include "Writer.hpp"

#include "Simulation/IfsfdMemory.hpp"
#include "Simulation/SfsfdMemory.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace VtkOutput {
template <typename TMemory>
std::unique_ptr<IterationResultWriter>
create_instance(TMemory const*                 memory) {
  return std::unique_ptr<IterationResultWriter>(new Writer<TMemory>(memory));
}

template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 0, 0, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 0, 1, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 1, 0, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 1, 1, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 0, 0, double, 3>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 0, 1, double, 3>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 1, 0, double, 3>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename SfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 1, 1, double, 3>
                ::MemoryType const*);

template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 0, 0, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 0, 1, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 1, 0, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 2>, 1, 1, double, 2>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 0, 0, double, 3>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 0, 1, double, 3>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 1, 0, double, 3>
                ::MemoryType const*);
template std::unique_ptr<IterationResultWriter>
create_instance(typename IfsfdSolverTraits
                <UniformGridGeometry<double, 3>, 1, 1, double, 3>
                ::MemoryType const*);

}
}
}
