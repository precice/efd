#pragma once

#include <boost/filesystem.hpp>

#include <string>

namespace FsiSimulation {
namespace FluidSimulation {
class IterationResultWriter {
public:
  virtual ~IterationResultWriter() {}

  virtual void
  setDestination(boost::filesystem::path const& directory_path,
                 std::string                    file_name_prefix) = 0;

  virtual void
  writeGeometry() = 0;

  virtual void
  writeAttributes() = 0;
};
}
}
