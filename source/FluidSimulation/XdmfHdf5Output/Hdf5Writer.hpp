#ifndef FsiSimulation_FluidSimulation_XdmfHdf5Output_Hdf5Writer_hpp
#define FsiSimulation_FluidSimulation_XdmfHdf5Output_Hdf5Writer_hpp

#include "FluidSimulation/BasicCell.hpp"
#include "FluidSimulation/ParallelDistribution.hpp"
#include "XdmfWriter.hpp"

#include <Eigen/Core>

#include <hdf5.h>

#include <boost/filesystem.hpp>

#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace FsiSimulation {
namespace FluidSimulation {
namespace XdmfHdf5Output {
template <typename TGrid>
class Hdf5Writer {
public:
  using GridType = TGrid;

  using CellAccessorType = typename GridType::CellAccessor;

  using CellType = typename CellAccessorType::CellType;

  using Scalar = typename CellType::Scalar;

  enum {
    Dimensions = CellType::Dimensions
  };

  using Path = boost::filesystem::path;

  using VectorDsize = Eigen::Matrix<hsize_t, Dimensions, 1>;

  using ParallelDistributionType = ParallelDistribution<Dimensions>;

  Hdf5Writer() {}

  void
  initialize(ParallelDistributionType const* parallel_distribution,
             GridType const*                 grid) {
    _parallelDistribution = parallel_distribution;
    _grid                 = grid;
    _corner               =
      _parallelDistribution->corner.template cast<hsize_t>();
    _globalSize =
      _parallelDistribution->globalCellSize.template cast<hsize_t>();
    _localSize =
      _parallelDistribution->localCellSize.template cast<hsize_t>();
  }

  VectorDsize const&
  size() const {
    return _globalSize;
  }

private:
  void
  initializeFile(std::string const& file_path,
                 hid_t&             file_id,
                 hid_t&             file_space,
                 hid_t&             chunk_space) const {
    MPI_Info info = MPI_INFO_NULL;
    // Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,
                     _parallelDistribution->mpiCommunicator,
                     info);
    // Create a new file collectively and release property list identifier.
    file_id = H5Fcreate(file_path.c_str(),
                        H5F_ACC_TRUNC,
                        H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    // Create the dataspace for the dataset.
    file_space = H5Screate_simple(Dimensions,
                                  _globalSize.data(),
                                  NULL);
    chunk_space = H5Screate_simple(Dimensions,
                                   _localSize.data(),
                                   NULL);
  }

private:
  void
  writeAttribute(hid_t const&       file_id,
                 hid_t const&       file_space,
                 hid_t const&       chunk_space,
                 std::string const& attribute_name,
                 hsize_t const&     size,
                 void*              data) {
    // Create chunked dataset.
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id,
                 Dimensions,
                 _localSize.data());
    hid_t dataset_id = H5Dcreate(file_id,
                                 (std::string("/") + attribute_name).c_str(),
                                 H5T_NATIVE_FLOAT,
                                 file_space,
                                 H5P_DEFAULT,
                                 plist_id,
                                 H5P_DEFAULT);
    H5Pclose(plist_id);

    auto size_vector = VectorDsize::Constant(size).eval();

    // Select hyperslab in the file.
    herr_t status = H5Sselect_hyperslab(file_space,
                                        H5S_SELECT_SET,
                                        _corner.data(),
                                        size_vector.data(),
                                        size_vector.data(),
                                        _localSize.data());

    // Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(dataset_id,
                      H5T_NATIVE_FLOAT,
                      chunk_space,
                      file_space,
                      plist_id,
                      data);

    // Close/release resources.
    H5Dclose(dataset_id);
    H5Pclose(plist_id);
  }

public:
  Path
  writeGeometry(Path const&        directory_path,
                std::string const& file_name_prefix,
                std::string const& dimensions_name) {
    Path current_file_path = directory_path;
    current_file_path.append(file_name_prefix + ".h5");

    hid_t file_id;
    hid_t file_space;
    hid_t chunk_space;

    initializeFile(current_file_path.string(),
                   file_id,
                   file_space,
                   chunk_space);

    float*      data       = new float[_localSize.prod() * Dimensions];
    std::size_t data_index = 0;

    for (auto const& accessor : _grid->innerGrid) {
      auto position = accessor.currentPosition().template cast<float>();

      for (int d = 0; d < Dimensions; ++d) {
        data[data_index] = position(d);
        ++data_index;
      }
    }

    writeAttribute(file_id,
                   file_space,
                   chunk_space,
                   dimensions_name,
                   Dimensions,
                   data);

    delete data;

    // Close/release resources.
    H5Sclose(chunk_space);
    H5Sclose(file_space);
    H5Fclose(file_id);

    return current_file_path;
  }

  Path
  writeAttributes(Path const&            directory_path,
                  std::string const& file_name_prefix,
                  std::vector<Attribute> attributes,
                  int const&             iteration_number) {
    Path current_file_path = directory_path;
    current_file_path.append(file_name_prefix
                             + "."
                             + std::to_string(iteration_number)
                             + ".h5");

    hid_t file_id;
    hid_t file_space;
    hid_t chunk_space;

    initializeFile(current_file_path.string(),
                   file_id,
                   file_space,
                   chunk_space);

    for (auto const& attribute : attributes) {
      int attribute_size = 1;

      if (attribute.type == Attribute::Type::Vector) {
        attribute_size = Dimensions;
      }

      float*      data       = new float[_localSize.prod()];
      std::size_t data_index = 0;

      for (auto const& accessor : _grid->innerGrid) {
        for (int j = 0; j < attribute_size; ++j) {
          data[data_index]
            = static_cast<float>(accessor.currentCell()->attribute(
                                   attribute.index, j));
        }
        ++data_index;
      }

      writeAttribute(file_id,
                     file_space,
                     chunk_space,
                     attribute.name,
                     attribute_size,
                     data);

      delete data;
    }

    H5Sclose(chunk_space);
    H5Sclose(file_space);
    H5Fclose(file_id);

    return current_file_path;
  }

private:
  ParallelDistributionType const* _parallelDistribution;
  GridType const*                 _grid;
  VectorDsize                     _corner;
  VectorDsize                     _globalSize;
  VectorDsize                     _localSize;
  Path                            _filePath;
};
}
}
}
#endif
