#include "Configuration.hpp"

#include <Uni/Firewall/Interface>
#include <Uni/Logging/format>

#include <boost/regex.hpp>

#include <map>

using Structure::Configuration;

namespace Structure {
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

Configuration::
Configuration() : _im(new ConfigurationImplementation(this)) {
  reset();
}

Configuration::
~Configuration() {}

void
Configuration::
reset() {

  _im->properties.clear();

  _set("/EnvironmentForce", false);
  _set("/Dimensions", false);
  _set("/PreciceConfigurationPath", false);
  _set("/Type", 0u);
  _set("/PreciceMode", false);
}

std::string
Configuration::
toString() const {
  // auto str = Uni::Logging::format(
  //   "Re = {1}\n"
  //    "Walls[2][1].velocity = {19}\n",
  //   any_to_string<std::string, long double>(
  //     _find("/Equations/Ins/ReynoldsNumber")),
  //   any_to_string<std::string, long double>(
  //     _find("/Equations/Ins/DiffusionMultiplier")),
  //   any_to_string<std::string, long double>(
  //     _find("/Equations/Ins/PressureGradientMultiplier")));

  return "";
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
