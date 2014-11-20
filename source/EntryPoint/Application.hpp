#ifndef FsiSimulation_EntryPoint_Arguments_hpp
#define FsiSimulation_EntryPoint_Arguments_hpp

#include <Uni/Firewall/Implementation>

namespace FsiSimulation {
namespace EntryPoint {
class ApplicationPrivateImplementation;
/**
 * Simple facility to drive the application
 */
class Application {
private:
  typedef ApplicationPrivateImplementation Implementation;
public:
  Application(int& argc, char** argv);

  Application(Application const& other) = delete;

  virtual
  ~Application();

  Application&
  operator=(Application const& other) = delete;

  void
  initialize();

  void
  run();

  void
  release();

private:
  void
  initializeOnMaster();

  void
  initializeOnSlave();

  void
  initializeOnMasterAndSlave();

  void
  runOnMaster();

  void
  runOnSlave();

  void
  runOnMasterAndSlave();

  void
  releaseOnMaster();

  void
  releaseOnSlave();

  void
  releaseOnMasterAndSlave();

  Uni_Firewall_IMPLEMENTATION_LINK(ApplicationPrivateImplementation)
};
}
}

#endif
