#ifndef FsiSimulation_GhostCellsHandler_hpp
#define FsiSimulation_GhostCellsHandler_hpp

#include "ParallelTopology.hpp"

#include <petscdm.h>
#include <petscdmda.h>

#include <array>
#include <functional>

namespace FsiSimulation {
template <typename TGrid, int D>
class GhostCellsHandler {
public:
  typedef ParallelTopology<D> SpecializedParallelTopology;
  typedef std::function<void (TGrid const*,
                              SpecializedParallelTopology const*,
                              Mat&)> PressureStencilHandler;
  template <typename THandler>
  using HandlerStack = std::array<std::array<THandler, 2>, D>;

  typedef HandlerStack<PressureStencilHandler> PressureStencilStack;

public:
  GhostCellsHandler() {}

  GhostCellsHandler(GhostCellsHandler const& other) = delete;

  ~GhostCellsHandler() {}

  GhostCellsHandler&
  operator=(GhostCellsHandler const& other) = delete;

  PressureStencilStack _pressureStencilStack;

private:
};
}
#endif
