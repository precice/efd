#ifndef FsiSimulation_FluidSimulation_VtkPlot_hpp
#define FsiSimulation_FluidSimulation_VtkPlot_hpp

#include "Grid.hpp"
#include "functions.hpp"

#include "ParallelDistribution.hpp"
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/locale.hpp>

#include <limits>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename TMemory,
          typename TGridGeometry,
          typename TScalar,
          int TD>
class VtkPlot {
public:
  typedef Grid<TMemory, TGridGeometry, TD>    GridType;
  typedef typename GridType::CellAccessorType CellAccessorType;
  typedef ParallelDistribution<TD>
    ParallelDistributionType;
  typedef typename CellAccessorType::GridGeometryType GridGeometry;

  typedef boost::filesystem::path Path;
  typedef boost::filesystem::fstream FileStream;
  typedef boost::locale::format      Format;

public:
  VtkPlot() {}

  void
  initialize(GridType const*                 grid,
             ParallelDistributionType const* parallelTopology,
             GridGeometry const*             gridGeometry,
             Path const&                     outputDirectory,
             std::string const&              fileNamePrefix) {
    _grid                 = grid;
    _parallelDistribution = parallelTopology;
    _gridGeometry         = gridGeometry;
    _outputDirectory      = outputDirectory;
    _fileNamePrefix       = fileNamePrefix;
    namespace bl          = boost::locale;
    _locale               = bl::generator().generate("en_US.UTF-8");
    _fileNamePrefix      += (Format(".{1}") % _parallelDistribution->rank)
                            .str(_locale);
  }

  void
  plot2(int const&     iterationCount,
        TScalar const& timeStamp,
        TScalar const& dt) {
    FileStream       fileStream;

    auto tempFileNamePrefix = _fileNamePrefix;
    tempFileNamePrefix += (Format(".{1}.my.vtk")
                           % iterationCount
                           // % timeStamp
                           ).str(_locale);

    auto tempPath = _outputDirectory / tempFileNamePrefix;

    fileStream.open(tempPath);

    if (!fileStream.is_open()) {
      fileStream.clear();
      fileStream.open(tempPath, FileStream::out);
    }
    // fileStream.close();

    fileStream << "# vtk DataFile Version 2.0" << std::endl
               << "I need something to put here" << std::endl
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
      pointsStream << accessor.currentPosition() (0) << " "
                   << accessor.currentPosition() (1);

      if (TD == 3) {
        pointsStream << " " << accessor.currentPosition() (2);
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
        rhsStream
          << PpeRhsGenerator::get(accessor, dt) << std::endl;
      } else {
        // rhsStream << "RHS "
        rhsStream
          << 0.0 << std::endl;
      }
      // fghStream << "FGH "
      fghStream
        << accessor.currentCell()->fgh(0) << " "
        << accessor.currentCell()->fgh(1);
      // velocityStream << "Velocity "
      velocityStream
        << accessor.currentCell()->velocity(0) << " "
        << accessor.currentCell()->velocity(1);

      if (TD == 3) {
        fghStream << " " << accessor.currentCell()->fgh(2);
        velocityStream << " " << accessor.currentCell()->velocity(2);
      }
      fghStream << std::endl;
      velocityStream << std::endl;
      // pressureStream << "Pressure "
      pressureStream
        << accessor.currentCell()->pressure() << std::endl;
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
  plot(int const&     iterationCount,
       TScalar const& timeStamp,
       TScalar const& dt) {
    FileStream       fileStream;

    auto tempFileNamePrefix = _fileNamePrefix;
    tempFileNamePrefix += (Format(".{1}.vtk")
                           % iterationCount
                           // % timeStamp
                           ).str(_locale);

    auto tempPath = _outputDirectory / tempFileNamePrefix;

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
                   % (TD == 2 ? " 1" : "")
                   % grid.innerSize().prod()).str(_locale);
    std::stringstream pointsStream;

    for (auto const& accessor : grid) {
      pointsStream << accessor.currentPosition() (0) << " "
                   << accessor.currentPosition() (1);

      if (TD == 3) {
        pointsStream << " " << accessor.currentPosition() (2);
      } else if (TD == 2) {
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
      velocityStream << accessor.currentCell()->velocity(0) << " "
                     << accessor.currentCell()->velocity(1);

      if (TD == 3) {
        velocityStream << " " << accessor.currentCell()->velocity(2);
      } else if (TD == 2) {
        velocityStream << " " << 0.0;
      }
      velocityStream << std::endl;
    }
    fileStream << velocityStream.str() << std::endl;

    fileStream << "\nSCALARS pressure float 1" << std::endl
               << "LOOKUP_TABLE default" << std::endl;

    for (auto const& accessor : _grid->innerGrid) {
      pressureStream << accessor.currentCell()->pressure() << std::endl;
    }
    fileStream << pressureStream.str() << std::endl;

    fileStream.close();
  }

private:
  GridType const*                 _grid;
  ParallelDistributionType const* _parallelDistribution;
  GridGeometry const*             _gridGeometry;
  Path                            _outputDirectory;
  std::locale                     _locale;
  std::string                     _fileNamePrefix;
};
}
}

#endif
