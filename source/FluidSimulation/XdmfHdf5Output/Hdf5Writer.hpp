#ifndef FsiSimulation_FluidSimulation_XdmfHdf5Output_Hdf5Writer_hpp
#define FsiSimulation_FluidSimulation_XdmfHdf5Output_Hdf5Writer_hpp

#include "FluidSimulation/BasicCell.hpp"
#include "FluidSimulation/ParallelDistribution.hpp"
#include "XdmfWriter.hpp"

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
  createFile(std::string const& file_path,
             hid_t&             file_id) const {
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
                 // Dimensions,
                 1,
                 &block_size);
    logInfo("{1}", attribute_name);
    hid_t dataset_id = H5Dcreate(file_id,
                                 (std::string("/") + attribute_name).c_str(),
                                 H5T_NATIVE_FLOAT,
                                 file_space,
                                 H5P_DEFAULT,
                                 plist_id,
                                 H5P_DEFAULT);
    H5Pclose(plist_id);

    logInfo("{1}", "sodfjslfj@@!!");

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

  void
  writeAttribute(hid_t const&       file_id,
                 std::string const& attribute_name,
                 hsize_t const&     data_size,
                 void*              data) {
    hsize_t size       = _globalSize.prod() * data_size;
    hsize_t start      = _corner.prod() * data_size;
    hsize_t count      = 1;
    hsize_t block_size = _localSize.prod() * data_size;

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
    logInfo("{1}", attribute_name);
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

    hid_t file_id;

    createFile(current_file_path.string(), file_id);

    for (int d = 0; d < dimension_names.size(); ++d) {
      float*  data = 0;
      hsize_t size = 1;

      if (d < Dimensions) {
        size = _localSize(d);

        if (_parallelDistribution->neighbors[d][1] < 0) {
          size += 1;
        }
        data = new float[size];
        std::size_t data_index = 0;

        auto accessor = *_grid->innerGrid.begin();

        for (int i = 0; i < size; ++i) {
          auto position = static_cast<float>(accessor.currentPosition()(d));
          ++accessor.index(d);
          data[data_index] = position;
          ++data_index;
        }
      } else {
        if (_parallelDistribution->rank != 0) {
          continue;
        }
        data    = new float[1];
        data[0] = 0.0;
      }
      // for (auto const& accessor : _grid->innerGrid) {
      // auto position = accessor.currentPosition().template cast<float>();

      // data[data_index] = position(d);
      // ++data_index;
      // }

      std::string name;

      writeGeometry(file_id,
                    dimension_names[d],
                    d,
                    size,
                    data);

      delete data;
    }

    // Close/release resources.
    H5Fclose(file_id);

    return current_file_path;
  }

  Path
  writeAttributes(Path const&            directory_path,
                  std::string const&     file_name_prefix,
                  std::vector<Attribute> attributes,
                  int const&             iteration_number) {
    Path current_file_path = directory_path;
    current_file_path.append(file_name_prefix
                             + "."
                             + std::to_string(iteration_number)
                             + ".h5");

    hid_t file_id;

    createFile(current_file_path.string(), file_id);

    for (auto const& attribute : attributes) {
      int attribute_size = 1;

      if (attribute.type == Attribute::Type::Vector) {
        attribute_size = 3;
      }

      float*      data       = new float[attribute_size * _localSize.prod()];
      std::size_t data_index = 0;

      for (auto const& accessor : _grid->innerGrid) {
        for (int j = 0; j < attribute_size; ++j) {
          if (j < Dimensions) {
            data[data_index]
              = static_cast<float>(accessor.currentCell()->attribute(attribute.index, j));
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
    }

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
