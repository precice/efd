#ifndef FsiSimulation_EntryPoint_VtkPlot_hpp
#define FsiSimulation_EntryPoint_VtkPlot_hpp

#include "Grid.hpp"
#include "ParallelTopology.hpp"

#include "stencils/mystencils.hpp"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/locale.hpp>

namespace FsiSimulation {
namespace EntryPoint {
template <typename TCellAccessor,
          typename Scalar,
          int D>
class VtkPlot {
public:
  typedef Grid<TCellAccessor, D>               SpecializedGrid;
  typedef ParallelTopology<D>                  SpecializedParallelTopology;
  typedef typename TCellAccessor::GridGeometry GridGeometry;

  typedef boost::filesystem::path Path;
  typedef
    boost::iostreams::stream<boost::iostreams::mapped_file_sink>
    MappedFileStream;
  typedef boost::filesystem::fstream FileStream;
  typedef boost::locale::format      Format;

public:
  VtkPlot() {}

  void
  initialize(SpecializedGrid const*             grid,
             SpecializedParallelTopology const* parallelTopology,
             GridGeometry const*                gridGeometry,
             Path const&                        outputDirectory,
             std::string const&                 fileNamePrefix) {
    _grid             = grid;
    _parallelTopology = parallelTopology;
    _gridGeometry     = gridGeometry;
    _outputDirectory  = outputDirectory;
    _fileNamePrefix   = fileNamePrefix;
    namespace bl      = boost::locale;
    _locale           = bl::generator().generate("en_US.UTF-8");
    _fileNamePrefix  += (Format(".{1}") % _parallelTopology->rank)
                        .str(_locale);
  }

  void
  plot2(int const&    iterationCount,
        Scalar const& timeStamp,
        Scalar const& dt) {
    MappedFileStream mappedFileStream;
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

    // mappedFileStream.open(tempPath.string());

    fileStream << "# vtk DataFile Version 2.0" << std::endl
               << "I need something to put here" << std::endl
               << "ASCII" << std::endl << std::endl;

    typedef typename SpecializedGrid::Base TempGrid;

    TempGrid grid = _grid->innerGrid;

    grid.setIndents(grid.leftIndent() - TempGrid::VectorDi::Ones(),
                    grid.rightIndent());

    fileStream << (Format("DATASET STRUCTURED_GRID\n"
                          "DIMENSIONS {1} 1\n"
                          "POINTS {2} float\n") %
                   grid.innerSize().transpose() %
                   grid.innerSize().prod()).str(_locale);
    std::stringstream pointsStream;

    for (auto const& accessor : grid) {
      // fileStream << (Format("{1}\n") %
      // accessor.currentPosition().transpose())
      // .str(_locale);
      pointsStream << accessor.currentPosition() (0) << " "
                   << accessor.currentPosition() (1) << " "
                   << 0 << std::endl;
      // logInfo("{1} {2}", accessor.indexValues().transpose(),
      // accessor.currentPosition().transpose());
    }
    fileStream << pointsStream.str() << std::endl;

    // fileStream << "\nCELL_DATA " << _grid->innerGrid.innerSize().prod()
    // << std::endl;

    // fileStream << "\nVECTORS velocity float" << std::endl;

    // for (auto const& accessor : _grid->innerGrid) {
    // fileStream << (Format("{1}\n")
    // % accessor.currentCell()->velocity().transpose())
    // .str(_locale);
    // }

    // fileStream << "\nVECTORS fgh float" << std::endl;

    std::stringstream rhsStream; // ! Stream for the pressure data
    std::stringstream fghStream; // ! Stream for the pressure data
    std::stringstream pressureStream; // ! Stream for the pressure data
    std::stringstream velocityStream; // ! Stream for the velocity data

    for (auto const& accessor :* _grid) {
      // typedef RhsProcessing
      // <typename SpecializedGrid::CellAccessor, Scalar, D> rhspr;
      // rhsStream << rhspr::compute(accessor, dt) << std::endl;
      fghStream << accessor.currentCell()->fgh(0) << " "
                << accessor.currentCell()->fgh(1) << " "
                << 0 << std::endl;
      velocityStream << accessor.currentCell()->velocity(0) << " "
                     << accessor.currentCell()->velocity(1) << " "
                     << 0 << std::endl;
      pressureStream << accessor.currentCell()->pressure() << std::endl;
    }

    fileStream << "CELL_DATA " << _grid->innerGrid.innerSize().prod()
               << std::endl
               << "SCALARS pressure float 1" << std::endl
               << "LOOKUP_TABLE default" << std::endl;
    fileStream << fghStream.str() << std::endl;
    // fileStream << rhsStream.str() << std::endl;
    fileStream << velocityStream.str() << std::endl;
    fileStream << pressureStream.str() << std::endl;

    // fileStream << "\nVECTORS velocity float" << std::endl;

    // for (auto const& accessor : _grid->innerGrid) {
    // fileStream << (Format("{1} {2}\n")
    // % accessor.indexValues().transpose()
    // % accessor.currentCell()->velocity().transpose())
    // .str(_locale);
    // }

    // fileStream << "\nSCALARS pressure float 1" << std::endl
    // << "LOOKUP_TABLE default" << std::endl;

    // for (auto const& accessor : _grid->innerGrid) {
    // fileStream << (Format("{1}\n") %
    // accessor.currentCell()->pressure()).str(_locale);
    // }

    fileStream.close();
  }

  void
  plot(int const&    iterationCount,
       Scalar const& timeStamp,
       Scalar const& dt) {
    MappedFileStream mappedFileStream;
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

    typedef typename SpecializedGrid::Base TempGrid;

    TempGrid grid = _grid->innerGrid;

    grid.setIndents(grid.leftIndent() - TempGrid::VectorDi::Ones(),
                    grid.rightIndent());

    fileStream << (Format("DATASET STRUCTURED_GRID\n"
                          "DIMENSIONS {1} 1\n"
                          "POINTS {2} float\n") %
                   grid.innerSize().transpose() %
                   grid.innerSize().prod()).str(_locale);
    std::stringstream pointsStream;

    for (auto const& accessor : grid) {
      pointsStream << accessor.currentPosition()(0) << " "
                   << accessor.currentPosition()(1) << " "
                   << "0.0" << std::endl;
    }
    fileStream << pointsStream.str() << std::endl;

    std::stringstream pressureStream;
    std::stringstream velocityStream;

    fileStream << "\nCELL_DATA " << _grid->innerGrid.innerSize().prod()
               << std::endl;

    fileStream << "\nVECTORS velocity float" << std::endl;

    for (auto const& accessor : _grid->innerGrid) {
      velocityStream << accessor.currentCell()->velocity(0) << " "
                     << accessor.currentCell()->velocity(1) << " "
                     << "0.0" << std::endl;
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
  SpecializedGrid const*             _grid;
  SpecializedParallelTopology const* _parallelTopology;
  GridGeometry const*                _gridGeometry;
  Path                               _outputDirectory;
  std::locale                        _locale;
  std::string                        _fileNamePrefix;
};
}
}

#endif
