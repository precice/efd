#pragma once

#include "FluidSimulation/functions.hpp"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/locale.hpp>

#include <limits>

namespace FsiSimulation {
namespace FluidSimulation {
namespace VtkOutput {
template <typename TMemory>
class Writer {
public:
  using MemoryType = TMemory;

  enum {
    Dimensions = MemoryType::Dimensions
  };

  using GridType = typename MemoryType::GridType;

  using Path = boost::filesystem::path;

  typedef boost::filesystem::fstream FileStream;
  typedef boost::locale::format      Format;

public:
  Writer() {}

  void
  initialize(MemoryType const* memory,
             Path const&       directory_path,
             std::string       file_name_prefix) {
    _memory         = memory;
    _directoryPath  = directory_path;
    _fileNamePrefix = file_name_prefix;
    _locale         = boost::locale::generator().generate("en_US.UTF-8");

    _fileNamePrefix
      += (Format(".{1}") % _memory->parallelDistribution()->rank).str(_locale);
  }

  void
  writeGeometry() {}

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

    std::stringstream pressureStream;
    std::stringstream velocityStream;

    fileStream << "\nCELL_DATA " <<
      _memory->grid()->innerGrid.innerSize().prod()
               << std::endl;

    fileStream << "\nVECTORS velocity float" << std::endl;

    for (auto const& accessor : _memory->grid()->innerGrid) {
      velocityStream << accessor.velocity(0) << " "
                     << accessor.velocity(1);

      if (Dimensions == 3) {
        velocityStream << " " << accessor.velocity(2);
      } else if (Dimensions == 2) {
        velocityStream << " " << 0.0;
      }
      velocityStream << std::endl;
    }
    fileStream << velocityStream.str() << std::endl;

    fileStream << "\nSCALARS pressure float 1" << std::endl
               << "LOOKUP_TABLE default" << std::endl;

    for (auto const& accessor : _memory->grid()->innerGrid) {
      pressureStream << accessor.pressure() << std::endl;
    }
    fileStream << pressureStream.str() << std::endl;

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
