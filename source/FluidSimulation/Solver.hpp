#pragma once

#include <boost/filesystem.hpp>

namespace precice {
class SolverInterface;
}

namespace FsiSimulation {
namespace FluidSimulation {
class Solver {
public:
  typedef boost::filesystem::path Path;

public:
  virtual
  ~Solver() {}

  virtual void
  initialize(precice::SolverInterface* preciceInteface,
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix) = 0;

  virtual bool
  iterate() = 0;
};
}
}
