#pragma once

#include <Uni/Firewall/Implementation>

namespace Fluid {
namespace EntryPoint {
class ApplicationPrivateImplementation;

/**
 * Simple facility to drive the application
 */
class Application {
private:
  using Implementation = ApplicationPrivateImplementation;

public:
  Application(int const& argc, char** argv);

  Application(Application const& other) = delete;

  virtual
  ~Application();

  Application&
  operator=(Application const& other) = delete;

  int
  parseArguments();

  void
  initialize();

  void
  run();

  void
  release();

private:
  void
  parseSimulationConfiguration();

  void
  createOutputDirectory();

  void
  initializePrecice();

  Uni_Firewall_IMPLEMENTATION_LINK(ApplicationPrivateImplementation)
};
}
}

