#include "Application.hpp"

#include "FluidSimulation/SimulationController.hpp"
#include "FluidSimulation/facade.hpp"
#include "XmlConfigurationParser.hpp"

#include "Utility/executablepath.hpp"
#include "Utility/pathoperations.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/Firewall/Interface>
#include <Uni/Logging/macros>

#include <boost/program_options.hpp>

#include "petscsys.h"

//

using FsiSimulation::EntryPoint::Application;

class FsiSimulation::EntryPoint::ApplicationPrivateImplementation {
  typedef precice::SolverInterface PreciceInterface;

  typedef std::unique_ptr<PreciceInterface> UniquePreciceInterface;

  typedef std::unique_ptr<FluidSimulation::SimulationController>
    UniqueSimulationController;

  typedef boost::filesystem::path Path;

  ApplicationPrivateImplementation(
    Application* in,
    int&         argc_,
    char**       argv_) : _in(in),
    argc(argc_),
    argv(argv_),
    masterRank(0) {
    namespace fs          = boost::filesystem;
    using LocaleGenerator = boost::locale::generator;
    globalLocale          = LocaleGenerator().generate("");
    std::locale::global(globalLocale);
    Path::imbue(std::locale());

    LocaleGenerator ansiLocaleGenerator;
    ansiLocaleGenerator.use_ansi_encoding(true);
    ansiLocale = ansiLocaleGenerator.generate("");

    applicationPath = Path(Utility::getExecutablePath()).parent_path();

    outputDirectoryPath    = fs::current_path();
    petscConfigurationPath = applicationPath;
    petscConfigurationPath.append("FluidPetsc/Basic.conf");
    preciceConfigurationPath = applicationPath;
    preciceConfigurationPath.append("Precice/SketchOfGeometryModeInFluid.xml");
    simulationConfigurationPath = applicationPath;
    simulationConfigurationPath.append("FluidSimulation/Channel.xml");
  }

  Uni_Firewall_INTERFACE_LINK(Application)

  int argc;
  char**      argv;
  std::locale globalLocale;
  std::locale ansiLocale;
  Path        applicationPath;

  Path        outputDirectoryPath;
  Path        vtkOutputDirectoryPath;
  Path        preciceConfigurationPath;
  Path        simulationConfigurationPath;
  Path        petscConfigurationPath;
  std::string vtkFilePrefix;

  int masterRank;
  int rank;
  int processCount;

  std::unique_ptr<FluidSimulation::Configuration> parameters;

  UniquePreciceInterface     preciceInterface;
  UniqueSimulationController mySimulation;

  bool
  isMaster() {
    return rank == masterRank;
  }
};

Application::
Application(int&   argc,
            char** argv)
  : _im(new Implementation(this,
                           argc,
                           argv)) {}

Application::
~Application() {}

int
Application::
parseArguments() {
  namespace po = boost::program_options;

  // Declare the supported options.
  po::options_description optionDescription("Allowed options");
  optionDescription.add_options()
    ("help,h", "Produce help message")
    ("output-directory,o",
    po::value<std::string>(),
    "Set output directory path")
    ("precice,p",
    po::value<std::string>(),
    "Set PreCiCe configuration path")
    ("simulation,s",
    po::value<std::string>(),
    "Set fluid simulation configuration path")
    ("petsc,e",
    po::value<std::string>(),
    "Set PETSc configuration path");

  po::variables_map options;
  po::store(po::parse_command_line(_im->argc, _im->argv, optionDescription),
            options);
  po::notify(options);

  if (options.count("help")) {
    std::cout << optionDescription << "\n";

    return 1;
  }

  if (options.count("output-directory")) {
    _im->outputDirectoryPath = options["output-directory"].as<std::string>();
    boost::filesystem::create_directories(_im->outputDirectoryPath);
    _im->outputDirectoryPath = boost::filesystem::canonical(
      options["output-directory"].as<std::string>());
  }

  if (options.count("petsc")) {
    _im->petscConfigurationPath = boost::filesystem::canonical(
      options["petsc"].as<std::string>());
  }

  if (options.count("precice")) {
    _im->preciceConfigurationPath = boost::filesystem::canonical(
      options["precice"].as<std::string>());
  }

  if (options.count("simulation")) {
    _im->simulationConfigurationPath = boost::filesystem::canonical(
      options["simulation"].as<std::string>());
  }

  return 0;
}

void
Application::
initialize() {
  // logInfo("Application path:\n{1}",
  // _im->applicationPath.string());
  // logInfo("Output directory path:\n{1}",
  // _im->outputDirectoryPath.string());
  // logInfo("PETSc configuration path:\n{1}",
  // _im->petscConfigurationPath.string());
  // logInfo("PreCICE configuration path:\n{1}",
  // _im->preciceConfigurationPath.string());
  // logInfo("Simulation configuration path:\n{1}",
  // _im->simulationConfigurationPath.string());

  auto petscConfigurationPath
    = Utility::convertUtfPathToAnsi(
    boost::filesystem::make_relative(
      _im->petscConfigurationPath).string(),
    _im->globalLocale,
    _im->ansiLocale);

  PetscInitialize(&_im->argc, &_im->argv,
                  petscConfigurationPath.c_str(),
                  PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &_im->processCount);
  MPI_Comm_rank(PETSC_COMM_WORLD, &_im->rank);

  parseSimulationConfiguration();
  createOutputDirectory();

  // Change current working directory of the application to overcome a Precice
  // configuration issue with python modules paths.
  boost::filesystem::current_path(_im->applicationPath);

  initializePrecice();

  _im->mySimulation
    = FluidSimulation::create_simulation_controller(_im->parameters.get());

  _im->mySimulation->initialize(_im->preciceInterface.get(),
                                _im->vtkOutputDirectoryPath,
                                _im->vtkFilePrefix);
}

void
Application::
run() {
  while (_im->mySimulation->iterate()) {
    if (_im->rank == 0) {
      logInfo("N = {1}; t = {2}",
              _im->mySimulation->iterationNumber(),
              _im->mySimulation->time());
    }
  }
}

void
Application::
release() {
  PetscFinalize();
}

void
Application::
parseSimulationConfiguration() {
  _im->parameters
    = XmlConfigurationParser::parse(_im->simulationConfigurationPath);
}

void
Application::
createOutputDirectory() {
  Implementation::Path outputDirectoryPath(_im->parameters->filename);

  if (outputDirectoryPath.is_relative()) {
    outputDirectoryPath = (_im->outputDirectoryPath / outputDirectoryPath);
  }

  auto outputFileName = outputDirectoryPath.filename();
  outputDirectoryPath = outputDirectoryPath.parent_path();
  boost::filesystem::create_directories(outputDirectoryPath);
  _im->vtkOutputDirectoryPath = outputDirectoryPath;
  outputDirectoryPath
    = boost::filesystem::make_relative(
    outputDirectoryPath) / outputFileName;
  _im->parameters->filename = outputDirectoryPath.string();
  _im->vtkFilePrefix        = outputFileName.string();
}

void
Application::
initializePrecice() {
  using namespace precice;

  _im->preciceInterface = Implementation::UniquePreciceInterface(
    new Implementation::PreciceInterface("Fluid",
                                         _im->rank,
                                         _im->processCount)
    );

  auto preciceConfigurationPath
    = Utility::convertUtfPathToAnsi(
    boost::filesystem::make_relative(
      _im->preciceConfigurationPath).string(),
    _im->globalLocale,
    _im->ansiLocale);

  _im->preciceInterface->configure(preciceConfigurationPath);
}
