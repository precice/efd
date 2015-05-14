#pragma once

#include "Simulation/IterationResultWriter.hpp"

#include <Uni/StructuredGrid/Basic/GlobalMultiIndex>
#include <Uni/StructuredGrid/Basic/Grid>

#include <boost/filesystem/fstream.hpp>
#include <boost/locale.hpp>

#include <limits>

namespace FsiSimulation {
namespace FluidSimulation {
namespace VtkOutput {
template <typename TMemory>
class Writer : public IterationResultWriter {
public:
  using MemoryType = TMemory;

  enum {
    Dimensions = MemoryType::Dimensions
  };

  using GridType = typename MemoryType::GridType;

  using CellAccessorType = typename MemoryType::CellAccessorType;

  using VectorDiType = typename MemoryType::VectorDiType;

  using ScalarType = typename MemoryType::ScalarType;

  using Path = boost::filesystem::path;

  typedef boost::filesystem::fstream FileStream;
  typedef boost::locale::format      Format;

  using AttributeType = typename MemoryType::AttributeType;

public:
  Writer(MemoryType const* memory) :
    _memory(memory) {}

  void
  setDestination(Path const& directory_path,
                 std::string file_name_prefix) {
    _directoryPath  = directory_path;
    _fileNamePrefix = file_name_prefix;
    _locale         = boost::locale::generator().generate("en_US.UTF-8");

    _fileNamePrefix
      += (Format(".{1}") % _memory->parallelDistribution()->rank()).str(_locale);
  }

  void
  writeAttributes() {
    FileStream fileStream;

    auto tempFileNamePrefix = _fileNamePrefix;
    tempFileNamePrefix += (Format(".{1}.vtk")
                           % _memory->iterationNumber()
                           ).str(_locale);

    auto tempPath = _directoryPath / tempFileNamePrefix;

    fileStream.open(tempPath);

    if (!fileStream.is_open()) {
      fileStream.clear();
      fileStream.open(tempPath, FileStream::out);
    }

    fileStream << "# vtk DataFile Version 2.0" << std::endl
               << "I need something to put here" << std::endl
               << "ASCII" << std::endl << std::endl;

    typedef typename GridType::SubgridType TempGrid;

    TempGrid point_grid = _memory->grid()->innerGrid;

    point_grid.setIndents(point_grid.leftIndent(),
                          point_grid.rightIndent() - TempGrid::VectorDi::Ones());

    fileStream << (Format("DATASET STRUCTURED_GRID\n"
                          "DIMENSIONS {1}{2}\n"
                          "POINTS {3} float\n")
                   % point_grid.innerSize().transpose()
                   % (Dimensions == 2 ? " 1" : "")
                   % point_grid.innerSize().prod()).str(_locale);
    std::stringstream pointsStream;

    for (auto const& accessor : point_grid) {
      pointsStream << accessor.position() (0) << " "
                   << accessor.position() (1);

      if (Dimensions == 3) {
        pointsStream << " " << accessor.position() (2);
      } else if (Dimensions == 2) {
        pointsStream << " " << 0.0;
      }
      pointsStream << std::endl;
    }
    fileStream << pointsStream.str() << std::endl;

    fileStream << "\nPOINT_DATA "
               << point_grid.innerSize().prod()
               << std::endl;

    int index = 0;

    for (auto const& attribute : * _memory->attributes()) {
      if (attribute.mode != AttributeType::DisplayMode::Points) {
        ++index;
        continue;
      }
      std::stringstream string_stream;
      struct GridTraits {
        using CellAccessorType
                = Uni::StructuredGrid::Basic::GlobalMultiIndex
                  <void, void, Dimensions, GridTraits>;

        using GridType
                = Uni::StructuredGrid::Basic::Grid
                  <void, void, Dimensions, GridTraits>;
      };

      std::array<typename GridTraits::GridType, Dimensions + 1> neighbors;
      std::array<ScalarType, Dimensions + 1>                    coeff;

      for (unsigned d = 0; d < Dimensions; ++d) {
        neighbors[d].initialize(VectorDiType::Constant(2),
                                VectorDiType::Zero(),
                                VectorDiType::Ones());
        coeff[d] = 1.0 / std::pow(2.0, Dimensions - 1);
      }
      neighbors[Dimensions].initialize(VectorDiType::Constant(2));
      coeff[Dimensions] = 1.0 / std::pow(2.0, Dimensions);

      if (attribute.type == AttributeType::Type::Vector) {
        fileStream << "\nVECTORS "
                   << attribute.name
                   << " float" << std::endl;

        for (auto const& accessor : point_grid) {
          for (int d = 0; d < 3; ++d) {
            if (d != 0) {
              string_stream << " ";
            }

            float temp_value = 0.0f;

            if (d < Dimensions) {
              ScalarType value = 0;

              for (auto const& accessor2 : neighbors[d]) {
                value += accessor.attribute(
                  accessor2.index() - VectorDiType::Ones(),
                  index,
                  d);
              }
              temp_value = static_cast<float>(coeff[d] * value);
            }
            string_stream << temp_value;
          }
          string_stream << std::endl;
        }
      } else {
        string_stream << "\nSCALARS "
                      << attribute.name
                      << " float 1" << std::endl
                      << "LOOKUP_TABLE default" << std::endl;

        for (auto const& accessor : point_grid) {
          ScalarType value = 0;

          for (auto const& accessor2 : neighbors[Dimensions]) {
            value += accessor.attribute(
              accessor2.index() - VectorDiType::Ones(),
              index);
          }
          string_stream << static_cast<float>(coeff[Dimensions] * value);
          string_stream << std::endl;
        }
      }
      fileStream << string_stream.str() << std::endl;
      ++index;
    }

    fileStream << "\nCELL_DATA "
               << _memory->grid()->innerGrid.innerSize().prod()
               << std::endl;

    index = 0;

    for (auto const& attribute : * _memory->attributes()) {
      if (attribute.mode != AttributeType::DisplayMode::Cells) {
        ++index;
        continue;
      }
      std::stringstream string_stream;

      if (attribute.type == AttributeType::Type::Vector) {
        fileStream << "\nVECTORS "
                   << attribute.name
                   << " float" << std::endl;

        for (auto const& accessor : _memory->grid()->innerGrid) {
          for (int d = 0; d < 3; ++d) {
            if (d != 0) {
              string_stream << " ";
            }

            float temp_value = 0.0f;

            if (d < Dimensions) {
              temp_value
                = static_cast<float>(accessor.centralizedAttribute(index, d));
            }
            string_stream << temp_value;
          }
          string_stream << std::endl;
        }
      } else {
        string_stream << "\nSCALARS "
                      << attribute.name
                      << " float 1" << std::endl
                      << "LOOKUP_TABLE default" << std::endl;

        for (auto const& accessor : _memory->grid()->innerGrid) {
          string_stream
            << static_cast<float>(accessor.centralizedAttribute(index));
          string_stream << std::endl;
        }
      }
      fileStream << string_stream.str() << std::endl;
      ++index;
    }

    fileStream.close();
  }

private:
  void
  computeCentredPointData(CellAccessorType const& accessor,
                          unsigned const&         attribute_index,
                          unsigned const&         dimension) const {
    VectorDiType offset = VectorDiType::Zero();
    ScalarType   value  = accessor.attribute(attribute_index, dimension);

    struct GridTraits {
      using CellAccessorType
              = Uni::StructuredGrid::Basic::GlobalMultiIndex
                <void, void, Dimensions, GridTraits>;

      using GridType
              = Uni::StructuredGrid::Basic::Grid
                <void, void, Dimensions, GridTraits>;
    };

    typename GridTraits::GridType neighbors;
    neighbors.initialize(VectorDiType::Constant(3));

    for (unsigned d = 0; d < Dimensions; ++d) {
      for (unsigned d2 = 0; d2 <= d; ++d2) {
        //
      }
    }
  }

  MemoryType const* _memory;
  std::locale       _locale;
  Path              _directoryPath;
  std::string       _fileNamePrefix;
};
}
}
}
