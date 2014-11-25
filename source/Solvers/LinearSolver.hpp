#ifndef FsiSimulation_Solvers_LenearSolver_hpp
#define FsiSimulation_Solvers_LenearSolver_hpp

#include "FlowField.h"
#include "Parameters.h"

namespace FsiSimulation {
namespace Solvers {
// Abstract class for linear solvers for the pressure
class LinearSolver {
protected:
  FlowField&        _flowField;  // ! Reference to the flow field
  const Parameters& _parameters; // ! Reference to the parameters

public:
  /** Constructor
   * @param flowField The flow field to solve
   * @param parameters Parameters of the problem
   */
  LinearSolver(FlowField& flowField, const Parameters& parameters);

  /** Solve the linear system for the pressure
   */
  virtual void
  solve() = 0;
};
}
}

#endif
