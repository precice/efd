#pragma once

#include <Uni/ExecutionControl/exception>
#include <Uni/Firewall/Implementation>

#include <Eigen/Core>

#include <boost/any.hpp>

#include <array>
#include <memory>
#include <string>

namespace Structure {
class ConfigurationImplementation;

class Configuration {
public:
  typedef long double                     ScalarType;

  typedef Eigen::Matrix<int, 3, 1>        VectorDiType;
  typedef Eigen::Matrix<ScalarType, 3, 1> VectorDsType;

public:
  Configuration();

  ~Configuration();

  void
  reset();

  std::string
  toString() const;

  bool
  is(std::string const& name) const;

  template <typename T>
  T
  get(std::string const& name) const {
    try {
      return boost::any_cast<T>(_find(name));
    } catch (boost::bad_any_cast const&) {
      logFatal("Cannot retrieve property '{1}' value", name);
      throw std::exception();
    }
  }

  template <typename T>
  T
  tryToGet(std::string const& name) const {
    return boost::any_cast<T>(_find(name));
  }

  template <typename T>
  bool
  isOfType(std::string const& name) const {
    auto value = _find(name);

    try {
      boost::any_cast<T>(value);

      return true;
    } catch (boost::bad_any_cast const&) {
      return false;
    }
  }

  void
  set(std::string const& name, boost::any const& value);

  void
  set(std::string const& name);

private:
  void
  _set(std::string const& name, boost::any const& value);

  boost::any const&
  _find(std::string const& name) const;

  boost::any&
  _find(std::string const& name);

  Uni_Firewall_IMPLEMENTATION_LINK(ConfigurationImplementation);
};
}
