#pragma once

#include "FluidSimulation/Configuration.hpp"

#include <boost/filesystem.hpp>

#include <memory>
#include <string>

namespace FsiSimulation {
namespace EntryPoint {
class XmlConfigurationParser {
public:
  static void
  parse(std::unique_ptr<FluidSimulation::Configuration> const& configuration,
        boost::filesystem::path const&                         filePath);

private:
  XmlConfigurationParser() {}

  XmlConfigurationParser(XmlConfigurationParser const&) = delete;

  ~XmlConfigurationParser() {}

  XmlConfigurationParser&
  operator=(XmlConfigurationParser const&) = delete;
};
}
}
