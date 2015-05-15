#include "Configuration.hpp"

#include <Uni/Firewall/Interface>
#include <Uni/Logging/format>

#include <boost/regex.hpp>

#include <map>

using FsiSimulation::FluidSimulation::Configuration;

namespace FsiSimulation {
namespace FluidSimulation {
template <typename T1, typename T2>
inline static std::string
any_to_string(boost::any const& any) {
  std::stringstream ss;
  try {
    auto temp = boost::any_cast<T1>(any);
    ss << temp;
  } catch (boost::bad_any_cast const&) {
    try {
      auto temp = boost::any_cast<T2>(any);
      ss << temp;
    } catch (boost::bad_any_cast const&) {
      throwException("Failed to cast any to string");
    }
  }

  return ss.str();
}

class ConfigurationImplementation {
  ConfigurationImplementation(Configuration* in) : _in(in) {}

  std::map<std::string, boost::any> properties;

  Uni_Firewall_INTERFACE_LINK(Configuration);
};
}
}

Configuration::
Configuration() : _im(new ConfigurationImplementation(this)) {
  reset();
}

Configuration::
~Configuration() {}

void
Configuration::
reset() {
  timeLimit           = 0;
  plotInterval        = 0;
  tau                 = 0.5;
  gamma               = 0.5;
  dimensions          = 2;
  width               = VectorDsType::Zero();
  size                = VectorDiType::Ones();
  parallelizationSize = VectorDiType::Ones();
  environment         = VectorDsType::Zero();
  solverType          = SolverEnum::Sfsfd;
  scalarType          = ScalarEnum::Double;
  outputType          = OutputEnum::Vtk;
  outerLayerSize      = 1;
  innerLayerSize      = 0;
  doImmersedBoundary  = false;
  doDebug             = false;

  for (int d = 0; d < 3; ++d) {
    walls[d][0].reset(new Input());
    walls[d][1].reset(new Input());
  }

  _im->properties.clear();

  _set("/Equations/Ins/ReynoldsNumber",
       (long double)(1));
  _set("/Equations/Ins/DiffusionMultiplier",
       std::string("default"));
  _set("/Equations/Ins/PressureGradientMultiplier",
       std::string("default"));

  _set("/Ib/PreciceConfigurationPath",           false);

  _set("/Ib/Features/DevelopingStructure",       false);

  _set("/Ib/Features/FullVelocityPrediction",    false);

  _set("/Ib/Schemes/DirectForcing/PreciceBased", false);
  _set(
    "/Ib/Schemes/DirectForcing/PreciceBased/OuterLayerSize",
    (unsigned)(0));
  _set(
    "/Ib/Schemes/DirectForcing/PreciceBased/InnerLayerSize",
    (unsigned)(1));

  _set("/Ib/Schemes/DirectForcing/RbfBased",
       false);
  _set("/Ib/Schemes/DirectForcing/RbfBased/SupportRadius",
       std::string("default"));
  _set("/Ib/Schemes/DirectForcing/RbfBased/ImqShape",
       std::string("default"));
}

std::string
Configuration::
toString() const {
  auto str = Uni::Logging::format(
    "Re = {1}\n"
    "DiffusionMultiplier = {2}\n"
    "PressureGradientMultiplier = {3}\n"
    "TimeLimit = {4}\n"
    "IterationLimit = {5}\n"
    "PlotInterval = {6}\n"
    "Tau = {7}\n"
    "Gamma = {8}\n"
    "Dimensions = {9}\n"
    "Width = {10}\n"
    "Size = {11}\n"
    "ParallelizationSize = {12}\n"
    "Environment = {13}\n"
    "Walls[0][0].velocity = {14}\n"
    "Walls[0][1].velocity = {15}\n"
    "Walls[1][0].velocity = {16}\n"
    "Walls[1][1].velocity = {17}\n"
    "Walls[2][0].velocity = {18}\n"
    "Walls[2][1].velocity = {19}\n",
    any_to_string<std::string, long double>(
      _find("/Equations/Ins/ReynoldsNumber")),
    any_to_string<std::string, long double>(
      _find("/Equations/Ins/DiffusionMultiplier")),
    any_to_string<std::string, long double>(
      _find("/Equations/Ins/PressureGradientMultiplier")),
    timeLimit,
    iterationLimit,
    plotInterval,
    tau,
    gamma,
    dimensions,
    width.transpose(),
    size.transpose(),
    parallelizationSize.transpose(),
    environment.transpose(),
    walls[0][0]->velocity().transpose(),
    walls[0][1]->velocity().transpose(),
    walls[1][0]->velocity().transpose(),
    walls[1][1]->velocity().transpose(),
    walls[2][0]->velocity().transpose(),
    walls[2][1]->velocity().transpose());

  str += Uni::Logging::format(
    "Do full IB velocity prediction = {1}\n",
    is("/Ib/Features/DevelopingStructure")
    ? "True" : "False");

  str += Uni::Logging::format(
    "Do full IB velocity prediction = {1}\n",
    is("/Ib/Features/FullVelocityPrediction")
    ? "True" : "False");

  if (is("/Ib/Schemes/DirectForcing/PreciceBased")) {
    str += Uni::Logging::format(
      "PreciceBased ImmersedBoundary Method\n"
      "OuterLayerSize = {1}\n"
      "InnerLayerSize = {2}\n",
      get<unsigned>(
        "/Ib/Schemes/DirectForcing/PreciceBased/OuterLayerSize"),
      get<unsigned>(
        "/Ib/Schemes/DirectForcing/PreciceBased/InnerLayerSize"));
  }

  return str;
}

bool
Configuration::
is(std::string const& name) const {
  auto find_it = _im->properties.find(name);

  if (find_it == _im->properties.end()) {
    throwException("There is no propery '{1}' to test", name);
  }

  try {
    return boost::any_cast<bool>(find_it->second);
  } catch (boost::bad_any_cast const&) {
    throwException("Cannot retrieve boolean property '{1}'", name);

    return false;
  }
}

void
Configuration::
set(std::string const& name, boost::any const&  value) {
  auto find_it = _im->properties.find(name);

  if (find_it == _im->properties.end()) {
    throwException("There is no propery '{1}' to set", name);
  }

  find_it->second = value;
}

void
Configuration::
set(std::string const& name) {
  set(name, true);
}

boost::any const&
Configuration::
_find(std::string const& name) const {
  auto find_it = _im->properties.find(name);

  if (find_it == _im->properties.end()) {
    throwException("There is no propery '{1}' to return", name);
  }

  return find_it->second;
}

boost::any&
Configuration::
_find(std::string const& name) {
  auto find_it = _im->properties.find(name);

  if (find_it == _im->properties.end()) {
    throwException("There is no propery '{1}' to return", name);
  }

  return find_it->second;
}

void
Configuration::
_set(std::string const& name, boost::any const&  value) {
  _im->properties.emplace(std::make_pair(name, value));
}
