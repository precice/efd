#ifndef FsiSimulation_GhostCellsHandler_hpp
#define FsiSimulation_GhostCellsHandler_hpp

#include "Solvers/GhostInitializationHandlers.hpp"
#include "Solvers/GhostMpiExchangeHandlers.hpp"
#include "Solvers/GhostPressureStencilHanlers.hpp"

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
  Solvers::Ghost::MpiExchange::FunctorStack<D>     _mpiVelocityExchangeStack;
  Solvers::Ghost::MpiExchange::FunctorStack<D>     _mpiPressureExchangeStack;
  Solvers::Ghost::Initialization::FunctorStack<D>  _fghInitialization;
  Solvers::Ghost::Initialization::FunctorStack<D>  _velocityInitialization;

private:
};
}
#endif
