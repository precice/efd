#ifndef FsiSimulation_FluidSimulation_XdmfHdf5Output_XdmfHdf5Writer_hpp
#define FsiSimulation_FluidSimulation_XdmfHdf5Output_XdmfHdf5Writer_hpp

#include "FluidSimulation/BasicCell.hpp"
#include "FluidSimulation/ParallelDistribution.hpp"
#include "Hdf5Writer.hpp"
#include "XdmfWriter.hpp"

#include <boost/filesystem.hpp>

#include <memory>
#include <string>
#include <vector>

namespace FsiSimulation {
namespace FluidSimulation {
namespace XdmfHdf5Output {
template <typename TGrid>
class XdmfHdf5Writer {
public:
  using GridType = TGrid;

  using CellAccessorType = typename GridType::CellAccessor;

  using CellType = typename CellAccessorType::CellType;

  enum {
    Dimensions = CellType::Dimensions
  };

  using Hdf5WriterType = Hdf5Writer<GridType>;

  using XdmfWriterType = XdmfWriter<Hdf5WriterType::Dimensions>;

  using UniqueXdmfWriterType = std::unique_ptr<XdmfWriterType>;

  using ParallelDistributionType
          = typename Hdf5WriterType::ParallelDistributionType;

  using Path = boost::filesystem::path;

  XdmfHdf5Writer() {}

  void
  initialize(ParallelDistributionType const* parallel_distribution,
             GridType const*                 grid,
             Path const&                     directory_path,
             std::string                     file_name_prefix) {
    _filePath       = directory_path;
    _fileNamePrefix = file_name_prefix;

    _hdf5Writer.initialize(parallel_distribution,
                           grid);

    if (parallel_distribution->rank == 0) {
      _xdmfWriter.reset(new XdmfWriterType());
    }

    int const attributes_size = CellType::Traits::getAttributesSize();

    for (int i = 0; i < attributes_size; ++i) {
      auto const& attribute = CellType::Traits::getAttribute(i);
      _attributes.emplace_back(attribute);
    }
  }

  void
  writeGeometry() {
    std::string dimensions_name;

    if (Dimensions == 2) {
      dimensions_name = "XY";
    } else {
      dimensions_name = "XYZ";
    }

    auto geometry_file_path = _hdf5Writer.writeGeometry(_filePath,
                                                        _fileNamePrefix,
                                                        dimensions_name);

    if (_xdmfWriter.get()) {
      _xdmfWriter->initialize(geometry_file_path.filename().string(),
                              dimensions_name,
                              _attributes,
                              _hdf5Writer.size().template cast<int>());
    }
  }

  void
  writeAttributes(int const& iteration_number) {
    auto attribute_file_path = _hdf5Writer.writeAttributes(
      _filePath,
      _fileNamePrefix,
      _attributes,
      iteration_number);

    if (_xdmfWriter.get()) {
      _xdmfWriter->write(
        _filePath,
        _fileNamePrefix,
        attribute_file_path.filename().string(),
        iteration_number);
    }
  }

private:
  Hdf5WriterType         _hdf5Writer;
  UniqueXdmfWriterType   _xdmfWriter;
  Path                   _filePath;
  std::string            _fileNamePrefix;
  std::vector<Attribute> _attributes;
};
}
}
}

#endif
