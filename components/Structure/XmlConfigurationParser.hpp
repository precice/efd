#pragma once

#include "Simulation.hpp"

#include <boost/filesystem.hpp>

#include <memory>

namespace Structure {
class XmlConfigurationParser {
public:
  static void
  parse(std::unique_ptr<Simulation> const& configuration,
        boost::filesystem::path const&                         filePath);

private:
  XmlConfigurationParser() {}

  XmlConfigurationParser(XmlConfigurationParser const&) = delete;

  ~XmlConfigurationParser() {}

  XmlConfigurationParser&
  operator=(XmlConfigurationParser const&) = delete;
};
}
