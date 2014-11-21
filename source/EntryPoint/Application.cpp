#include "Application.hpp"

#include "Configuration.hpp"
#include "FlowField.h"
#include "MeshsizeFactory.hpp"
#include "Parameters.h"
#include "Simulation.h"
#include "TurbulentFlowField.h"
#include "TurbulentSimulation.h"
#include "parallelManagers/PetscParallelConfiguration.h"

#include <Uni/ExecutionControl/exception>
#include <Uni/Firewall/Interface>
#include <Uni/Logging/logging>
#include <Uni/Logging/logging>
#include <Uni/Platform/operatingsystem>

#include <boost/filesystem.hpp>
#include <boost/locale.hpp>

#include <mpi.h>
#include <petscsys.h>

#include <iostream>
#include <memory>

#if defined (Uni_Platform_OS_WINDOWS)
#  include <windows.h>
#elif defined (Uni_Platform_OS_UNIX)
#  include <unistd.h>
#else
#  error "Unknown platform"
#endif
//

using FsiSimulation::EntryPoint::Application;

namespace FsiSimulation {
namespace EntryPoint {
#if defined (Uni_Platform_OS_WINDOWS)
boost::filesystem::path
_getExecPath() {
  TCHAR buffer[MAX_PATH];
  auto  size =  GetModuleFileNameW(NULL, buffer, MAX_PATH);

  if (size == MAX_PATH ||
      size == 0) {
    throwException("Failed to locate application");
  }
  namespace fs = boost::filesystem;
  fs::path result(buffer);

  return result;
}
#elif defined (Uni_Platform_OS_UNIX)
boost::filesystem::path
_getExecPath() {
  char buffer[PATH_MAX];
  auto size = readlink("/proc/self/exe", buffer, PATH_MAX);

  if (size == -1) {
    throwException("Failed to locate application");
  }
  namespace fs = boost::filesystem;
  fs::path result(buffer);

  return result;
}
#else
#  error "Unknown platform"
#endif
}
}

class FsiSimulation::EntryPoint::ApplicationPrivateImplementation {
  typedef std::unique_ptr<Simulation> UniqueSimulation;
  typedef std::unique_ptr<FlowField>  UniqueFlowField;
  typedef boost::filesystem::path     Path;

  ApplicationPrivateImplementation(
    Application* in,
    int&         argc_,
    char**       argv_) : _in(in),
                          argc(argc_),
                          argv(argv_),
                          masterRank(0) {
    using LocaleGenerator = boost::locale::generator;
    globalLocale          = LocaleGenerator().generate("");
    // Create and install global locale
    std::locale::global(globalLocale);
    // Make boost.filesystem use it
    Path::imbue(std::locale());

    LocaleGenerator ansiLocaleGenerator;
    ansiLocaleGenerator.use_ansi_encoding(true);
    ansiLocale = ansiLocaleGenerator.generate("");

    applicationPath = _getExecPath();

    petscConfigurationPath = applicationPath.parent_path();
    petscConfigurationPath.append("FluidPetsc/Basic.conf");

    logInfo("Utf-8 locale: {1}",
            std::use_facet<boost::locale::info>(globalLocale).name());
    logInfo("Ansi locale: {1}",
            std::use_facet<boost::locale::info>(ansiLocale).name());
    logInfo("Application path:\n{1}",         applicationPath.string());
    logInfo("PETSc configuration path:\n{1}", petscConfigurationPath.string());
  }

  Uni_Firewall_INTERFACE_LINK(Application)

  int    argc;
  char**      argv;
  std::locale globalLocale;
  std::locale ansiLocale;
  Path        applicationPath;
  Path        petscConfigurationPath;

  int masterRank;
  int rank;
  int processCount;

  Parameters parameters;

  UniqueSimulation simulation;
  UniqueFlowField  flowField;

  bool
  isMaster() { return rank == masterRank; }
};

Application::
Application(int&   argc,
            char** argv)
  : _im(new Implementation(this,
                           argc,
                           argv)) {
  auto toEncoding = std::use_facet<boost::locale::info>(
    _im->ansiLocale).encoding();
  auto fromEncoding = std::use_facet<boost::locale::info>(
    _im->globalLocale).encoding();
  auto petscConfigurationPath =
    boost::locale::conv::between(_im->petscConfigurationPath.string(),
                                 toEncoding, fromEncoding);
  PetscInitialize(&argc, &argv,
                  petscConfigurationPath.c_str(),
                  PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &_im->processCount);
  MPI_Comm_rank(PETSC_COMM_WORLD, &_im->rank);
}

Application::
~Application() {}

void
Application::
initialize() {
  if (_im->isMaster()) {
    initializeOnMaster();
  } else {
    initializeOnSlave();
  }
}

void
Application::
run() {
  if (_im->isMaster()) {
    runOnMaster();
  } else {
    runOnSlave();
  }
}

void
Application::
release() {
  if (_im->isMaster()) {
    releaseOnMaster();
  } else {
    releaseOnSlave();
  }
}

void
Application::
initializeOnMaster() {
  initializeOnMasterAndSlave();
}

void
Application::
initializeOnSlave() {
  initializeOnMasterAndSlave();
}

void
Application::
initializeOnMasterAndSlave() {
  Configuration configuration(_im->argv[1]);
  configuration.loadParameters(_im->parameters);
  PetscParallelConfiguration parallelConfiguration(_im->parameters);
  MeshsizeFactory::getInstance().initMeshsize(_im->parameters);

  // #ifndef NDEBUG
  std::cout << "Processor " << _im->parameters.parallel.rank << " with index ";
  std::cout << _im->parameters.parallel.indices[0] << ",";
  std::cout << _im->parameters.parallel.indices[1] << ",";
  std::cout << _im->parameters.parallel.indices[2];
  std::cout <<    " is computing the size of its subdomain and obtains ";
  std::cout << _im->parameters.parallel.localSize[0] << ", ";
  std::cout << _im->parameters.parallel.localSize[1] << " and ";
  std::cout << _im->parameters.parallel.localSize[2];
  std::cout << ". Left neighbour: " << _im->parameters.parallel.leftNb;
  std::cout << ", Right neighbour: " << _im->parameters.parallel.rightNb;
  std::cout << ". Bottom neighbour: " << _im->parameters.parallel.bottomNb;
  std::cout << ", Top neighbour: " << _im->parameters.parallel.topNb;
  std::cout << ". Back neighbour: " << _im->parameters.parallel.backNb;
  std::cout << ", Front neighbour: " << _im->parameters.parallel.frontNb;
  std::cout << std::endl;
  std::cout << "Min. meshsizes: " << _im->parameters.meshsize->getDxMin() <<
    ", " << _im->parameters.meshsize->getDyMin() << ", " <<
    _im->parameters.meshsize->getDzMin() <<
    std::endl;
  // #endif

  // initialise simulation
  if (_im->parameters.simulation.type == "turbulence") {
    if (_im->isMaster()) {
      std::cout << "Start turbulence simulation in " <<
        _im->parameters.geometry.dim << "D" << std::endl;
    }
    auto turbulentFlowField = new TurbulentFlowField(_im->parameters);

    _im->flowField  = Implementation::UniqueFlowField(turbulentFlowField);
    _im->simulation =
      Implementation::UniqueSimulation(
        new TurbulentSimulation(_im->parameters, *turbulentFlowField));
  } else if (_im->parameters.simulation.type == "dns") {
    if (_im->isMaster()) {
      std::cout << "Start DNS simulation in " << _im->parameters.geometry.dim <<
        "D" << std::endl;
    }
    auto flowField = new FlowField(_im->parameters);

    _im->flowField  = Implementation::UniqueFlowField(flowField);
    _im->simulation =
      Implementation::UniqueSimulation(
        new Simulation(_im->parameters, *flowField));
  } else {
    handleError(1,
                "Unknown simulation type! Currently supported: dns, turbulence");
  }

  _im->simulation->initializeFlowField();
  // flowField->getFlags().show();
}

void
Application::
runOnMaster() {
  runOnMasterAndSlave();
}

void
Application::
runOnSlave() {
  runOnMasterAndSlave();
}

void
Application::
runOnMasterAndSlave() {
  FLOAT time       = 0.0;
  FLOAT timeVtk    = _im->parameters.vtk.interval;
  FLOAT timeStdOut = _im->parameters.stdOut.interval;
  int   timeSteps  = 0;

  // plot initial state
  _im->simulation->plotVTK(timeSteps);

  // time loop
  while (time < _im->parameters.simulation.finalTime) {
    _im->simulation->solveTimestep();

    time += _im->parameters.timestep.dt;

    // std-out: terminal info
    if (_im->isMaster() && (timeStdOut <= time)) {
      std::cout << "Current time: " << time << "\ttimestep: " <<
        _im->parameters.timestep.dt << std::endl;
      timeStdOut += _im->parameters.stdOut.interval;
    }

    // VTK output
    if (timeVtk <= time) {
      _im->simulation->plotVTK(timeSteps);
      timeVtk += _im->parameters.vtk.interval;
    }
    timeSteps++;
  }

  // plot final output
  _im->simulation->plotVTK(timeSteps);
}

void
Application::
releaseOnMaster() {
  releaseOnMasterAndSlave();
}

void
Application::
releaseOnSlave() {
  releaseOnMasterAndSlave();
}

void
Application::
releaseOnMasterAndSlave() {
  PetscFinalize();
}
