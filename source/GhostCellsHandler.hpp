#ifndef FsiSimulation_GhostCellsHandler_hpp
#define FsiSimulation_GhostCellsHandler_hpp

#include "Solvers/GhostPressureStencilHanlers.hpp"
#include "Solvers/GhostMpiExchangeHandlers.hpp"

namespace FsiSimulation {
template <int D>
class GhostCellsHandler {
public:
  GhostCellsHandler() {}

  GhostCellsHandler(GhostCellsHandler const& other) = delete;

  ~GhostCellsHandler() {}

  GhostCellsHandler&
  operator=(GhostCellsHandler const& other) = delete;

  Solvers::Ghost::PressureStencil::FunctorStack<D> _pressureStencilStack;
  Solvers::Ghost::MpiExchange::FunctorStack<D> _mpiVelocityExchangeStack;

private:
};
}
#endif
