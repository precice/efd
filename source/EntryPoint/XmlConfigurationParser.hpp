#pragma once

#include "Simulation/Configuration.hpp"

#include <Uni/Firewall/Implementation>

#include <boost/filesystem.hpp>

#include <memory>

namespace Fluid {
namespace EntryPoint {
class XmlConfigurationParserImplementation;
class XmlConfigurationParser {
public:
  XmlConfigurationParser(
    std::unique_ptr<Simulation::Configuration> const& configuration,
    boost::filesystem::path const&                         filePath);

  XmlConfigurationParser(XmlConfigurationParser const&) = delete;

  ~XmlConfigurationParser();

  XmlConfigurationParser&
  operator=(XmlConfigurationParser const&) = delete;

private:
  Uni_Firewall_IMPLEMENTATION_LINK(XmlConfigurationParserImplementation);
};
}
}
