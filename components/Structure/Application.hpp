#pragma once

#include <Uni/Firewall/Implementation>

namespace Structure {
class ApplicationPrivateImplementation;

class Application {
private:
  typedef ApplicationPrivateImplementation Implementation;

public:
  Application(int const& argc, char** argv);

  Application(Application const&) = delete;

  ~Application();

  Application&
  operator=(Application const&) = delete;

  int
  parseArguments();

  void
  initialize();

  void
  run();

private:
  void
  initializePrecice();

  Uni_Firewall_IMPLEMENTATION_LINK(ApplicationPrivateImplementation)
};
}
