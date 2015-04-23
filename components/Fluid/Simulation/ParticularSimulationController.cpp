#include "ParticularSimulationController.hpp"

#include "FsfdSolver.hpp"
#include "IfsfdMemory.hpp"
#include "IterationResultWriter.hpp"
#include "Reporter.hpp"
#include "SfsfdMemory.hpp"

#include <Uni/Firewall/Interface>
#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TSolverTraits>
class ParticularSimulationControllerImplementation {
public:
  using Interface = ParticularSimulationController<TSolverTraits>;

  ParticularSimulationControllerImplementation(Interface* in) : _in(in) {}

  typename Interface::SolverType                             solver;
  std::unique_ptr<IterationResultWriter> resultWriter;
  Reporter*                              reporter;

  long double lastTime;
  long double timeLimit;
  long double lastPlotTimeStamp;
  typename Interface::ScalarType         plotInterval;
  unsigned long long iterationLimit;

  Uni_Firewall_INTERFACE_LINK(ParticularSimulationController<TSolverTraits>);
};

template <typename T>
ParticularSimulationController<T>::
ParticularSimulationController() : _im(new Implementation(this)) {}

template <typename T>
ParticularSimulationController<T>::
~ParticularSimulationController() {}

template <typename T>
void
ParticularSimulationController<T>::
setIterationResultWriter(IterationResultWriter* iteration_result_writer) {
  _im->resultWriter.reset(iteration_result_writer);
}

template <typename T>
typename ParticularSimulationController<T>::SolverType const*
ParticularSimulationController<T>::
solver() const {
  return &_im->solver;
}

template <typename T>
typename ParticularSimulationController<T>::SolverType*
ParticularSimulationController<T>::
solver() {
  return &_im->solver;
}

template <typename T>
long double const&
ParticularSimulationController<T>::
time() const {
  return _im->solver.memory()->time();
}

template <typename T>
long double const&
ParticularSimulationController<T>::
timeLimit() const {
  return _im->timeLimit;
}

template <typename T>
long double&
ParticularSimulationController<T>::
timeLimit() {
  return _im->timeLimit;
}

template <typename T>
long double const&
ParticularSimulationController<T>::
lastPlotTimeStamp() const {
  return _im->lastPlotTimeStamp;
}

template <typename T>
long double&
ParticularSimulationController<T>::
lastPlotTimeStamp() {
  return _im->lastPlotTimeStamp;
}

template <typename T>
typename ParticularSimulationController<T>::ScalarType const&
ParticularSimulationController<T>::
plotInterval() const {
  return _im->plotInterval;
}

template <typename T>
typename ParticularSimulationController<T>::ScalarType&
ParticularSimulationController<T>::
plotInterval() {
  return _im->plotInterval;
}

template <typename T>
unsigned long long const&
ParticularSimulationController<T>::
iterationNumber() const {
  return _im->solver.memory()->iterationNumber();
}

template <typename T>
unsigned long long const&
ParticularSimulationController<T>::
iterationLimit() const {
  return _im->iterationLimit;
}

template <typename T>
unsigned long long&
ParticularSimulationController<T>::
iterationLimit() {
  return _im->iterationLimit;
}

template <typename T>
void
ParticularSimulationController<T>::
initialize(precice::SolverInterface* preciceInteface,
           Reporter*                 reporter,
           Path const&               outputDirectory,
           std::string const&        fileNamePrefix) {
  // logInfo("Enter Initialize");
  _im->reporter = reporter;

  _im->reporter->setAt("Name", fileNamePrefix);

  // logInfo("Result Writer is being set its destination");
  _im->resultWriter->setDestination(outputDirectory, fileNamePrefix);

  // logInfo("Solver is being initialized");
  _im->solver.initialize(preciceInteface, reporter);

  _im->resultWriter->initialize();

  _im->lastPlotTimeStamp = 0.0;
  _im->lastTime          = -1.0;

  if (_im->plotInterval >= 0) {
    _im->resultWriter->writeGeometry();
    _im->resultWriter->writeAttributes();
  }
}

template <typename T>
bool
ParticularSimulationController<T>::
iterate() {
  if (std::abs(_im->lastTime - _im->solver.memory()->time())
      <= std::numeric_limits<long double>::epsilon()) {
    return false;
  }

  if (_im->timeLimit > 0
      && _im->solver.memory()->time() >= _im->timeLimit) {
    return false;
  }

  if (_im->iterationLimit > 0
      && _im->solver.memory()->iterationNumber() >= _im->iterationLimit) {
    return false;
  }
  _im->lastTime = _im->solver.memory()->time();

  _im->solver.iterate();

  if (_im->plotInterval >= 0) {
    if ((_im->solver.memory()->time() - _im->lastPlotTimeStamp)
        > _im->plotInterval) {
      _im->lastPlotTimeStamp = _im->solver.memory()->time();

      _im->resultWriter->writeAttributes();
    }
  }

  return true;
}

Fluid_InstantiateExternTemplates(ParticularSimulationController);
}
}
