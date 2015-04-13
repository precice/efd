#pragma once

#include "Simulation/IterationResultWriter.hpp"

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
      += (Format(".{1}") % _memory->parallelDistribution()->rank).str(_locale);
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

    typedef typename GridType::BaseType TempGrid;

    TempGrid grid = _memory->grid()->innerGrid;

    grid.setIndents(grid.leftIndent(),
                    grid.rightIndent() - TempGrid::VectorDi::Ones());

    fileStream << (Format("DATASET STRUCTURED_GRID\n"
                          "DIMENSIONS {1}{2}\n"
                          "POINTS {3} float\n")
                   % grid.innerSize().transpose()
                   % (Dimensions == 2 ? " 1" : "")
                   % grid.innerSize().prod()).str(_locale);
    std::stringstream pointsStream;

    for (auto const& accessor : grid) {
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

    fileStream << "\nCELL_DATA "
               << _memory->grid()->innerGrid.innerSize().prod()
               << std::endl;

    int index = 0;

    for (auto const& attribute : * _memory->attributes()) {
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
  MemoryType const* _memory;
  std::locale       _locale;
  Path              _directoryPath;
  std::string       _fileNamePrefix;
};
}
}
}
