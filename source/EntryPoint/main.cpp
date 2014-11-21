#include "Configuration.hpp"
#include "MeshsizeFactory.hpp"
#include "Simulation.h"
#include "TurbulentSimulation.h"
#include "parallelManagers/PetscParallelConfiguration.h"

#include "Application.hpp"

#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

int
main(int argc, char* argv[]) {
  using namespace FsiSimulation::EntryPoint;

  Application application(argc, argv);

  auto result = application.parseArguments();
  if (result == 1) {
    return 0;
  }

  if (result != 0) {
    return result;
  }

  application.initialize();

  application.run();

  application.release();
}
