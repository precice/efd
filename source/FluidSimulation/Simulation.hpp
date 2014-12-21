#ifndef FsiSimulation_FluidSimulation_Simulation_hpp
#define FsiSimulation_FluidSimulation_Simulation_hpp

#include <precice/SolverInterface.hpp>

#include <boost/filesystem.hpp>

namespace FsiSimulation {
namespace FluidSimulation {
class Simulation {
public:
  typedef boost::filesystem::path Path;

public:
  virtual
  ~Simulation() {}

  virtual void
  initialize(precice::SolverInterface* preciceInteface,
             Path const&               outputDirectory,
             std::string const&        fileNamePrefix) = 0;

  virtual bool
  iterate() = 0;
};
}
}
#endif
