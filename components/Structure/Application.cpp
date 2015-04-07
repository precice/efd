#include "Application.hpp"

#include "Simulation.hpp"
#include "XmlConfigurationParser.hpp"

#include <precice/SolverInterface.hpp>

#include <Uni/ExecutionControl/exception>
#include <Uni/Firewall/Interface>
#include <Uni/Logging/macros>
#include <Uni/executablepath>
#include <Uni/pathoperations>

#include <boost/filesystem.hpp>
#include <boost/locale.hpp>
#include <boost/program_options.hpp>

#include <memory>

using Structure::Application;

class Structure::ApplicationPrivateImplementation {
  using PreciceInterface = precice::SolverInterface;

  using UniquePreciceInterface = std::unique_ptr<PreciceInterface>;

  using Path = boost::filesystem::path;

  ApplicationPrivateImplementation(Application* in,
                                   int const&   argc_,
                                   char**       argv_) : _in(in),
    argc(argc_),
    argv(argv_),
    simulation(new Simulation()) {
    namespace fs          = boost::filesystem;
    using LocaleGenerator = boost::locale::generator;
    globalLocale          = LocaleGenerator().generate("");
    std::locale::global(globalLocale);
    Path::imbue(std::locale());

    LocaleGenerator ansiLocaleGenerator;
    ansiLocaleGenerator.use_ansi_encoding(true);
    ansiLocale = ansiLocaleGenerator.generate("");

    applicationPath = Path(Uni::getExecutablePath()).parent_path();
  }

  Uni_Firewall_INTERFACE_LINK(Application)

  int argc;
  char**      argv;
  std::locale globalLocale;
  std::locale ansiLocale;
  Path        applicationPath;

  Path preciceConfigurationPath;
  Path simulationConfigurationPath;

  UniquePreciceInterface      preciceInterface;
  std::unique_ptr<Simulation> simulation;
};

Application::
Application(int const& argc, char** argv)
  : _im(new Implementation(this, argc, argv)) {}

Application::
~Application() {
  //
}

int
Application::
parseArguments() {
  namespace po = boost::program_options;

  po::options_description optionDescription("Allowed options");
  optionDescription.add_options()
    ("help,h", "Produce help message")
    ("precice,p",
    po::value<std::string>(),
    "Set PreCiCe configuration path")
    ("simulation,s",
    po::value<std::string>(),
    "Set simulation configuration path");

  po::variables_map options;
  po::store(po::parse_command_line(_im->argc, _im->argv, optionDescription),
            options);
  po::notify(options);

  if (options.count("help")) {
    std::cout << optionDescription << "\n";

    return 1;
  }

  if (options.count("precice")) {
    _im->preciceConfigurationPath
      = boost::filesystem::canonical(options["precice"].as<std::string>());
  } else {
    throwException("No precice configuration provided, try --help, -h");
  }

  if (options.count("simulation")) {
    _im->simulationConfigurationPath = options["simulation"].as<std::string>();
  } else {
    throwException("No configuration was provided, try --help, -h");
  }

  return 0;
}

void
Application::
initialize() {
  XmlConfigurationParser::parse(
    _im->simulation,
    _im->simulationConfigurationPath);

  // Change current working directory of the application to overcome a Precice
  // configuration issue with python modules paths.
  boost::filesystem::current_path(_im->preciceConfigurationPath.parent_path());

  initializePrecice();

  _im->simulation->initialize(_im->preciceInterface.get());
}

void
Application::
run() {
  while (_im->simulation->iterate()) {
    //
  }
  _im->preciceInterface->finalize();
}

void
Application::
initializePrecice() {
  using namespace precice;

  _im->preciceInterface.reset(
    new Implementation::PreciceInterface("Structure", 0, 1));

  auto preciceConfigurationPath
    = Uni::convertUtfPathToAnsi(
    boost::filesystem::make_relative(
      _im->preciceConfigurationPath).string(),
    _im->globalLocale,
    _im->ansiLocale);

  _im->preciceInterface->configure(preciceConfigurationPath);
}
