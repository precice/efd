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

  application.initialize();

  application.run();

  application.release();
}
