#pragma once

namespace FsiSimulation {
namespace FluidSimulation {
class LinearSolver {
public:
  virtual ~LinearSolver() {}

  virtual long double
  absoluteTolerance() const = 0;

  virtual void
  absoluteTolerance(long double const& tolerance) = 0;

  virtual long double
  relativeTolerance() const = 0;

  virtual void
  relativeTolerance(long double const& tolerance) = 0;

  virtual void
  update() = 0;

  virtual void
  solve() = 0;

  virtual void
  release() = 0;
};
}
}
