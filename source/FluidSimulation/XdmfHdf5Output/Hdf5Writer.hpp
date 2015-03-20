#pragma once

#include <Uni/Logging/macros>

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
template <typename TMemory>
class Hdf5Writer {
public:
  using MemoryType = TMemory;

  using GridType = typename MemoryType::GridType;

  using CellAccessorType = typename GridType::CellAccessor;

  enum {
    Dimensions = CellAccessorType::Dimensions
  };

  using ScalarType = typename MemoryType::ScalarType;

  using ParallelDistributionType
          = typename MemoryType::ParallelDistributionType;

  using AttributeType = typename MemoryType::AttributeType;

  using Path = boost::filesystem::path;

  using VectorDsize = Eigen::Matrix<hsize_t, Dimensions, 1>;

  Hdf5Writer() {}

  void
  initialize(MemoryType const* memory) {
    _memory = memory;
    _corner =
      _memory->parallelDistribution()->corner.template cast<hsize_t>();
    _globalSize =
      _memory->parallelDistribution()->globalCellSize.template cast<hsize_t>();
    _localSize =
      _memory->parallelDistribution()->localCellSize.template cast<hsize_t>();
  }

  VectorDsize const&
  size() const {
    return _globalSize;
  }

private:
  void
  createFile(std::string const& file_path,
             hid_t&             file_id) const {
    MPI_Info info = MPI_INFO_NULL;
    // Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,
                     _memory->parallelDistribution()->mpiCommunicator,
                     info);
    // Create a new file collectively and release property list identifier.
    file_id = H5Fcreate(file_path.c_str(),
                        H5F_ACC_TRUNC,
                        H5P_DEFAULT,
                        plist_id);
    H5Pclose(plist_id);
  }

private:
  void
  writeGeometry(hid_t const&       file_id,
                std::string const& attribute_name,
                int const&         dimension,
                hsize_t const&     block_size,
                void*              data) {
    hsize_t size  = block_size;
    hsize_t start = 0;
    hsize_t count = 1;

    if (dimension < Dimensions) {
      size  = _globalSize(dimension) + 1;
      start = _corner(dimension);
    }
    // Create the dataspace for the dataset.
    hid_t file_space = H5Screate_simple(
      1,
      &size,
      NULL);
    hid_t chunk_space = H5Screate_simple(
      1,
      &block_size,
      NULL);

    // Create chunked dataset.
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id,
                 1,
                 &block_size);
    hid_t dataset_id = H5Dcreate(file_id,
                                 (std::string("/") + attribute_name).c_str(),
                                 H5T_NATIVE_FLOAT,
                                 file_space,
                                 H5P_DEFAULT,
                                 plist_id,
                                 H5P_DEFAULT);
    H5Pclose(plist_id);

    // Select hyperslab in the file.
    herr_t status = H5Sselect_hyperslab(
      file_space,
      H5S_SELECT_SET,
      &start,
      &count,
      &count,
      &block_size);

    // Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    status = H5Dwrite(dataset_id,
                      H5T_NATIVE_FLOAT,
                      chunk_space,
                      file_space,
                      plist_id,
                      data);

    // Close/release resources.
    H5Dclose(dataset_id);
    H5Pclose(plist_id);
    H5Sclose(chunk_space);
    H5Sclose(file_space);
  }

  void
  writeAttribute(hid_t const&       file_id,
                 std::string const& attribute_name,
                 hsize_t const&     data_size,
                 void*              data) {
    VectorDsize size;
    VectorDsize start;
    VectorDsize count;
    VectorDsize block_size;

    for (int d = 0; d < Dimensions - 1; ++d) {
      size(d)       = _globalSize(Dimensions - 1 - d);
      start(d)      = _corner(Dimensions - 1 - d);
      count(d)      = 1;
      block_size(d) = _localSize(Dimensions - 1 - d);
    }

    size(Dimensions - 1)       = _globalSize(0) * data_size;
    start(Dimensions - 1)      = _corner(0) * data_size;
    count(Dimensions - 1)      = 1;
    block_size(Dimensions - 1) = _localSize(0) * data_size;

    // Create the dataspace for the dataset.
    hid_t file_space = H5Screate_simple(
      Dimensions,
      size.data(),
      NULL);
    hid_t chunk_space = H5Screate_simple(
      Dimensions,
      block_size.data(),
      NULL);

    hid_t plist_id;

    // Create chunked dataset.
    // hid_t plist_id = H4Pcreate(H5P_DATASET_CREATE);
    // H5Pset_chunk(plist_id,
    // Dimensions,
    // block_size.data());
    //
    hid_t dataset_id = H5Dcreate(file_id,
                                 (std::string("/") + attribute_name).c_str(),
                                 H5T_NATIVE_FLOAT,
                                 file_space,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);
    // H5Pclose(plist_id);

    // Select hyperslab in the file.
    herr_t status = H5Sselect_hyperslab(
      file_space,
      H5S_SELECT_SET,
      start.data(),
      count.data(),
      count.data(),
      block_size.data());

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
    H5Sclose(chunk_space);
    H5Sclose(file_space);
  }

public:
  Path
  writeGeometry(Path const&                     directory_path,
                std::string const&              file_name_prefix,
                std::vector<std::string> const& dimension_names) {
    Path current_file_path = directory_path;
    current_file_path.append(file_name_prefix + ".h5");

    if (_memory->parallelDistribution()->rank != 0) {
      return current_file_path;
    }

    hid_t file_id = H5Fcreate(current_file_path.c_str(),
                              H5F_ACC_TRUNC,
                              H5P_DEFAULT,
                              H5P_DEFAULT);

    for (int d = 0; d < dimension_names.size(); ++d) {
      float*  data = 0;
      hsize_t size = 0;

      if (d < Dimensions) {
        size = _globalSize(d) + 1;

        data = new float[size];
        std::size_t data_index = 0;

        auto accessor = *_memory->grid()->innerGrid.begin();

        for (int i = 0; i < size; ++i) {
          auto position = static_cast<float>(
            accessor.memory()->gridGeometry()->computeCellPosition(d, i));
          data[data_index] = position;
          ++data_index;
        }
      } else {
        size    = 1;
        data    = new float[1];
        data[0] = 0.0;
      }

      hid_t file_space = H5Screate_simple(1, &size, NULL);
      hid_t dataset_id = H5Dcreate(
        file_id,
        (std::string("/") + dimension_names[d]).c_str(),
        H5T_NATIVE_FLOAT,
        file_space,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT);
      herr_t status = H5Dwrite(dataset_id,
                               H5T_NATIVE_FLOAT,
                               H5S_ALL,
                               H5S_ALL,
                               H5P_DEFAULT,
                               data);

      delete data;
    }

    // Close/release resources.
    H5Fclose(file_id);

    return current_file_path;
  }

  Path
  writeAttributes(Path const&        directory_path,
                  std::string const& file_name_prefix) {
    Path current_file_path = directory_path;
    current_file_path.append(
      file_name_prefix
      + "."
      + std::to_string(_memory->iterationNumber())
      + ".h5");

    hid_t file_id;

    createFile(current_file_path.string(), file_id);

    int attribute_index = 0;

    for (auto const& attribute : * _memory->attributes()) {
      int attribute_size = 1;

      if (attribute.type == AttributeType::Type::Vector) {
        attribute_size = 3;
      }

      float*      data       = new float[attribute_size * _localSize.prod()];
      std::size_t data_index = 0;

      for (auto const& accessor : _memory->grid()->innerGrid) {
        for (int j = 0; j < attribute_size; ++j) {
          if (j < Dimensions) {
            data[data_index]
             = static_cast<float>(accessor.attribute(attribute_index, j));
          } else {
            data[data_index] = 0.0;
          }
          ++data_index;
        }
      }

      writeAttribute(file_id,
                     attribute.name,
                     attribute_size,
                     data);

      delete data;
      ++attribute_index;
    }

    H5Fclose(file_id);

    return current_file_path;
  }

private:
  MemoryType const* _memory;
  VectorDsize       _corner;
  VectorDsize       _globalSize;
  VectorDsize       _localSize;
  Path              _filePath;
};
}
}
}
