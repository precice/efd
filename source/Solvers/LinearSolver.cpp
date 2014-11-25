#include "LinearSolver.hpp"
//

using FsiSimulation::Solvers::LinearSolver;

LinearSolver::
LinearSolver(FlowField& flowField, const Parameters& parameters)
  : _flowField(flowField), _parameters(parameters) {}
