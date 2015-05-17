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
    structureConfiguration(new Configuration()) {
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
  Path structureConfigurationPath;

  std::unique_ptr<Configuration> structureConfiguration;
  UniquePreciceInterface         preciceInterface;
  std::unique_ptr<Simulation>    simulation;
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
    ("structure,s",
    po::value<std::string>(),
    "Set structure configuration path");

  po::variables_map options;
  po::store(po::parse_command_line(_im->argc, _im->argv, optionDescription),
            options);
  po::notify(options);

  if (options.count("help")) {
    std::cout << optionDescription << "\n";

    return 1;
  }

  if (options.count("structure")) {
    _im->structureConfigurationPath = options["structure"].as<std::string>();
  } else {
    throwException("No configuration was provided, try --help, -h");
  }

  return 0;
}

void
Application::
initialize() {
  XmlConfigurationParser(_im->structureConfiguration,
                         _im->structureConfigurationPath);

  initializePrecice();

  _im->simulation.reset(new Simulation(_im->structureConfiguration.get())),

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

  auto config_path
    = _im->structureConfiguration->get<boost::filesystem::path>(
    "/PreciceConfigurationPath");

  config_path = boost::filesystem::canonical(config_path);

  logInfo("{1}", config_path);

  // Change current working directory of the application to overcome a Precice
  // configuration issue with python modules paths.
  boost::filesystem::current_path(config_path.parent_path());

  _im->preciceInterface.reset(
    new Implementation::PreciceInterface("Structure", 0, 1));

  auto preciceConfigurationPath
    = Uni::convertUtfPathToAnsi(
    boost::filesystem::make_relative(config_path).string(),
    _im->globalLocale, _im->ansiLocale);

  _im->preciceInterface->configure(preciceConfigurationPath);
}
