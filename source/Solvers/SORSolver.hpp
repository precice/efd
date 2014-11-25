#ifndef FsiSimulation_Solvers_SORSolver_hpp
#define FsiSimulation_Solvers_SORSolver_hpp

// This moduke solves the diffusion equation for the pressure. Right now it's
// very simple.

#include "../Definitions.h"
#include "../FlowField.h"
#include "../Parameters.h"
#include "LinearSolver.hpp"

namespace FsiSimulation {
namespace Solvers {
class SORSolver : public LinearSolver {
public:
  SORSolver(FlowField& flowField, const Parameters& parameters);

  void
  solve();
};
}
}

#endif
