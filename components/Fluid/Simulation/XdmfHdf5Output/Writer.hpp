#pragma once

#include "Simulation/IterationResultWriter.hpp"

#include "Private/Hdf5Writer.hpp"
#include "Private/XdmfWriter.hpp"

#include <boost/filesystem.hpp>

#include <memory>
#include <string>
#include <vector>

namespace FsiSimulation {
namespace FluidSimulation {
namespace XdmfHdf5Output {
template <typename TMemory>
class Writer : public IterationResultWriter {
public:
  using MemoryType = TMemory;

  using Hdf5WriterType = Private::Hdf5Writer<MemoryType>;

  using XdmfWriterType = Private::XdmfWriter<MemoryType>;

  using UniqueXdmfWriterType = std::unique_ptr<XdmfWriterType>;

  using Path = boost::filesystem::path;

  Writer(MemoryType const* memory) : _memory(memory) {}

  void
  setDestination(Path const& directory_path,
                 std::string file_name_prefix) {
    _directoryPath  = directory_path;
    _fileNamePrefix = file_name_prefix;
  }

  void
  initialize() {
    _timeStepFileNames.clear();

    _hdf5Writer.initialize(_memory);

    if (_memory->parallelDistribution()->rank == 0) {
      _xdmfWriter.reset(new XdmfWriterType());
      _xdmfWriter->initialize(_memory);
    }
  }

  void
  writeGeometry() {
    std::vector<std::string> dimension_names;

    dimension_names.push_back("X");
    dimension_names.push_back("Y");
    dimension_names.push_back("Z");

    auto geometry_file_path = _hdf5Writer.writeGeometry(_directoryPath,
                                                        _fileNamePrefix,
                                                        dimension_names);

    if (_xdmfWriter.get()) {
      _xdmfWriter->writeGeometry(geometry_file_path.filename().string(),
                                 dimension_names,
                                 _hdf5Writer.size().template cast<int>());
    }
  }

  void
  writeAttributes() {
    auto attribute_file_path = _hdf5Writer.writeAttributes(_directoryPath,
                                                           _fileNamePrefix);

    if (_xdmfWriter.get()) {
      Path xmdf_file_path = _directoryPath;
      xmdf_file_path.append(_fileNamePrefix
                            + "."
                            + std::to_string(_memory->iterationNumber())
                            + ".xdmf");
      _timeStepFileNames.push_back(xmdf_file_path.filename().string());

      _xdmfWriter->writeAttributes(xmdf_file_path,
                                   attribute_file_path.filename().string());

      xmdf_file_path = _directoryPath;
      xmdf_file_path.append(_fileNamePrefix + ".xdmf");

      _xdmfWriter->writeTemporal(xmdf_file_path,
                                 _timeStepFileNames);
    }
  }

private:
  MemoryType const*        _memory;
  Hdf5WriterType           _hdf5Writer;
  UniqueXdmfWriterType     _xdmfWriter;
  Path                     _directoryPath;
  std::string              _fileNamePrefix;
  std::vector<std::string> _timeStepFileNames;
};
}
}
}
