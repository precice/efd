#ifndef FsiSimulation_EntryPoint_XmlConfigurationParser_hpp
#define FsiSimulation_EntryPoint_XmlConfigurationParser_hpp

#include "FluidSimulation/Configuration.hpp"

#include <boost/filesystem.hpp>

#include <memory>
#include <string>

namespace FsiSimulation {
namespace EntryPoint {
class XmlConfigurationParser {
public:
  static
  std::unique_ptr<FluidSimulation::Configuration>
  parse(boost::filesystem::path const& filePath);

private:
  XmlConfigurationParser() {}

  XmlConfigurationParser(XmlConfigurationParser const&) = delete;

  ~XmlConfigurationParser() {}

  XmlConfigurationParser&
  operator=(XmlConfigurationParser const&) = delete;
};
}
}

#endif
