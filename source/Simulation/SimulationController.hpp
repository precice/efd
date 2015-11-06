#pragma once

#include <boost/filesystem.hpp>

namespace precice {
class SolverInterface;
}

namespace Fluid {
namespace Simulation {
class Reporter;

class SimulationController {
public:
  typedef boost::filesystem::path Path;

public:
  virtual ~SimulationController() {}

  virtual long double const&
  time() const = 0;

  virtual unsigned long long const&
  iterationNumber() const = 0;

  virtual void
  initialize(precice::SolverInterface* preciceInteface,
             Reporter*                 reporter,
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix) = 0;

  virtual bool
  iterate() = 0;
};
}
}
