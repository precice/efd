#ifndef FsiSimulation_EntryPoint_Configuration_hpp
#define FsiSimulation_EntryPoint_Configuration_hpp

#include "FluidSimulation/Configuration.hpp"

#include <petscksp.h>
#include <string>

namespace FsiSimulation {
namespace EntryPoint {
class Configuration {
private:
  std::string _filename;
  int         _dim;

public:
  Configuration();

  Configuration(std::string const& filename);

  void
  setFileName(std::string const& filename);

  void
  loadParameters(FluidSimulation::Configuration& parameters,
                 MPI_Comm const&                 communicator =
                   PETSC_COMM_WORLD);
};
}
}

#endif
