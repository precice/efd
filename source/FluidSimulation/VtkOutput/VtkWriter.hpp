#ifndef FsiSimulation_FluidSimulation_VtkOutput_VtkWriter_hpp
#define FsiSimulation_FluidSimulation_VtkOutput_VtkWriter_hpp

#include "FluidSimulation/ParallelDistribution.hpp"
#include "FluidSimulation/functions.hpp"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/locale.hpp>

#include <limits>

namespace FsiSimulation {
namespace FluidSimulation {
namespace VtkOutput {
template <typename TGrid>
class VtkWriter {
public:
  using GridType = TGrid;

  using CellAccessorType = typename GridType::CellAccessorType;

  using CellType = typename CellAccessorType::CellType;

  using Scalar = typename CellType::Scalar;

  enum {
    Dimensions = CellType::Dimensions
  };

  typedef ParallelDistribution<Dimensions> ParallelDistributionType;

  typedef boost::filesystem::path    Path;
  typedef boost::filesystem::fstream FileStream;
  typedef boost::locale::format      Format;

public:
  VtkWriter() {}

  void
  initialize(ParallelDistributionType const* parallel_distribution,
             GridType const*                 grid,
             Path const&                     directory_path,
             std::string                     file_name_prefix) {
    _grid                 = grid;
    _parallelDistribution = parallel_distribution;
    _directoryPath        = directory_path;
    _fileNamePrefix       = file_name_prefix;
    _locale               = boost::locale::generator().generate("en_US.UTF-8");
    _fileNamePrefix      += (Format(".{1}") % _parallelDistribution->rank)
                            .str(_locale);
  }

  void
  writeGeometry() {}

  void
  plot2(int const&    iterationCount,
        Scalar const& timeStamp,
        Scalar const& dt) {
    FileStream fileStream;

    auto tempFileNamePrefix = _fileNamePrefix;
    tempFileNamePrefix += (Format(".{1}.my.vtk")
                           % iterationCount
                           // % timeStamp
                           ).str(_locale);

    auto tempPath = _directoryPath / tempFileNamePrefix;

    fileStream.open(tempPath);

    if (!fileStream.is_open()) {
      fileStream.clear();
      fileStream.open(tempPath, FileStream::out);
    }
    // fileStream.close();

    fileStream << "# vtk DataFile Version 2.0" << std::endl
               << "Fsi Simulation Results" << std::endl
               << "ASCII" << std::endl << std::endl;

    typedef typename GridType::Base TempGrid;

    TempGrid grid = _grid->innerGrid;

    grid.setIndents(grid.leftIndent(),
                    grid.rightIndent() - TempGrid::VectorDi::Ones());

    fileStream << (Format("DATASET STRUCTURED_GRID\n"
                          "DIMENSIONS {1} \n"
                          "POINTS {2} float\n") %
                   grid.innerSize().transpose() %
                   grid.innerSize().prod()).str(_locale);
    std::stringstream pointsStream;

    for (auto const& accessor : grid) {
      pointsStream << accessor.position() (0) << " "
                   << accessor.position() (1);

      if (Dimensions == 3) {
        pointsStream << " " << accessor.position() (2);
      }
      pointsStream << std::endl;
    }
    fileStream << pointsStream.str() << std::endl;

    std::stringstream rhsStream;
    std::stringstream fghStream;
    std::stringstream pressureStream;
    std::stringstream velocityStream;

    fghStream.precision(std::numeric_limits<float>::digits10);
    rhsStream.precision(std::numeric_limits<float>::digits10);
    velocityStream.precision(std::numeric_limits<float>::digits10);
    pressureStream.precision(std::numeric_limits<float>::digits10);

    for (auto const& accessor : _grid->innerGrid) {
      if (accessor.indexValues()
          != CellAccessorType::VectorDiType::Zero() &&
          (accessor.indexValues().array()
           <= (_grid->innerGrid.innerLimit().array())).all()) {
        // rhsStream << "RHS "
        // rhsStream
        // << PpeRhsGenerator::get(accessor, dt) << std::endl;
      } else {
        // rhsStream << "RHS "
        rhsStream
          << 0.0 << std::endl;
      }
      // fghStream << "FGH "
      fghStream
        << accessor.fgh(0) << " "
        << accessor.fgh(1);
      // velocityStream << "Velocity "
      velocityStream
        << accessor.velocity(0) << " "
        << accessor.velocity(1);

      if (Dimensions == 3) {
        fghStream << " " << accessor.fgh(2);
        velocityStream << " " << accessor.velocity(2);
      }
      fghStream << std::endl;
      velocityStream << std::endl;
      // pressureStream << "Pressure "
      pressureStream
        << accessor.pressure() << std::endl;
    }

    fileStream << "CELL_DATA " << _grid->innerGrid.innerSize().prod()
               << std::endl;
    fileStream << fghStream.str() << std::endl;
    fileStream << rhsStream.str() << std::endl;
    fileStream << velocityStream.str() << std::endl;
    fileStream << pressureStream.str() << std::endl;

    fileStream.close();
  }

  void
  writeAttributes(int const&    iteration_number,
                  double const& time) {
    ((void)time);
    FileStream fileStream;

    auto tempFileNamePrefix = _fileNamePrefix;
    tempFileNamePrefix += (Format(".{1}.vtk")
                           % iteration_number
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

    typedef typename GridType::Base TempGrid;

    TempGrid grid = _grid->innerGrid;

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

    fileStream << "\nCELL_DATA " << _grid->innerGrid.innerSize().prod()
               << std::endl;

    fileStream << "\nVECTORS velocity float" << std::endl;

    for (auto const& accessor : _grid->innerGrid) {
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

    for (auto const& accessor : _grid->innerGrid) {
      pressureStream << accessor.pressure() << std::endl;
    }
    fileStream << pressureStream.str() << std::endl;

    fileStream.close();
  }

  void
  simplePlot(int const&    iterationCount,
             Scalar const& timeStamp,
             Scalar const& dt) {
    FileStream fileStream;

    auto tempFileNamePrefix = _fileNamePrefix;
    tempFileNamePrefix += (Format(".{1}.vtk")
                           % iterationCount
                           // % timeStamp
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

    typedef typename GridType::Base TempGrid;

    TempGrid grid = _grid->innerGrid;

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

    fileStream << "\nCELL_DATA " << _grid->innerGrid.innerSize().prod()
               << std::endl;

    fileStream << "\nVECTORS velocity float" << std::endl;

    for (auto const& accessor : _grid->innerGrid) {
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

    for (auto const& accessor : _grid->innerGrid) {
      pressureStream << accessor.pressure() << std::endl;
    }
    fileStream << pressureStream.str() << std::endl;

    fileStream.close();
  }

private:
  GridType const*                 _grid;
  ParallelDistributionType const* _parallelDistribution;
  Path                            _directoryPath;
  std::locale                     _locale;
  std::string                     _fileNamePrefix;
};
}
}
}

#endif
