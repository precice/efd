#include "SfsfdHandlersBuilder.hpp"

namespace FsiSimulation {
namespace FluidSimulation {
namespace GhostLayer {
template <typename T, int U, int D>
SfsfdHandlersBuilder<T, U, D>::
SfsfdHandlersBuilder(Configuration* configuration, SolverType*    simulation)
  : BaseType(configuration, simulation) {}

template <typename T, int U, int D>
void
SfsfdHandlersBuilder<T, U, D>::
setAsInput() {
  this->BaseType::setAsInput();
}

template <typename T, int U, int D>
void
SfsfdHandlersBuilder<T, U, D>::
setAsParabolicInput() {
  this->BaseType::setAsParabolicInput();
}

template <typename T, int U, int D>
void
SfsfdHandlersBuilder<T, U, D>::
setAsOutput() {
  this->BaseType::setAsOutput();
}

template <typename T, int U, int D>
void
SfsfdHandlersBuilder<T, U, D>::
setAsMpiExchange() {
  this->BaseType::setAsMpiExchange();
}

template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 0, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 1, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 0, double, 2>, 2, 1>;

template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 0, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 1, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 0, 1, double, 2>, 2, 1>;

template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 0, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 1, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 0, double, 2>, 2, 1>;

template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 0, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 1, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 2>, 1, 1, double, 2>, 2, 1>;

template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 0, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 1, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 0, double, 3>, 2, 1>;

template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 0, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 1, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 0, 1, double, 3>, 2, 1>;

template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 0, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 1, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 0, double, 3>, 2, 1>;

template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 0, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 1, 1>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 0>;
template class SfsfdHandlersBuilder
  <SfsfdSolverTraits<UniformGridGeometry<double, 3>, 1, 1, double, 3>, 2, 1>;
}
}
}
