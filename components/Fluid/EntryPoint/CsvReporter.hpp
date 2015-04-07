#pragma once

#include "Simulation/Reporter.hpp"

#include <Uni/ExecutionControl/assert>

#include <boost/filesystem/fstream.hpp>
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include <iostream>
#include <vector>

namespace FsiSimulation {
namespace EntryPoint {
class CsvReporter : public FluidSimulation::Reporter {
private:
  using FileMapping = boost::interprocess::file_mapping;

  using MappedRegion = boost::interprocess::mapped_region;

  using FileStream = boost::filesystem::fstream;

public:
  using Path = boost::filesystem::path;

  CsvReporter() : _isInitialized(false) {}

  CsvReporter(CsvReporter const&) = delete;

  ~CsvReporter() {}

  CsvReporter&
  operator=(CsvReporter const&) = delete;

  void
  initialize(Path const& directory_path,
             std::string file_name_prefix) {
    release();

    _infoFilePath = directory_path;
    _infoFilePath.append(file_name_prefix + "-info.csv");

    _filePath = directory_path;
    _filePath.append(file_name_prefix + ".csv");

    _chunkSize = MappedRegion::get_page_size();

    std::size_t multiplier = 1;

    while (multiplier * _chunkSize < 1000) {
      ++multiplier;
    }
    _chunkSize = multiplier * _chunkSize;

    _fileSize     = 0;
    _position     = 0;
    _regionOffset = 0;

    resizeFileAndMapRegion();

    _isInitialized = true;
  }

  void
  setAt(unsigned const&    column,
        std::string const& column_name,
        std::string const& value) {
    if (_info.size() <= column) {
      _info.resize(column + 1);
    }

    _info[column] = std::make_pair(column_name, value);
  }

  void
  addAt(unsigned const& column, std::string const& value) {
    if (_row.size() <= column) {
      _row.resize(column + 1);
    }

    _row[column] = value;
  }

  void
  append(std::string const& value) {
    std::size_t new_position = _position + value.size();

    if (new_position >= _region.get_size()) {
      resizeFileAndMapRegion(value.size());
      new_position = _position;
    }

    char* data = static_cast<char*>(_region.get_address());
    data += _position;
    std::memcpy(data, value.data(), value.size());
    _position = new_position;
  }

  void
  recordInfo() {
    std::string result = createRow(_info);
    _info.clear();

    FileStream fileStream(_infoFilePath,
                          std::fstream::binary
                          | std::fstream::trunc
                          | std::fstream::out);
    fileStream << result;
    fileStream.close();
  }

  void
  recordIteration() {
    std::string result = createRow(_row);
    _row.clear();

    result += '\n';
    this->append(result);
  }

  void
  release() {
    if (_isInitialized) {
      boost::filesystem::resize_file(_filePath, _regionOffset + _position);
      _region.flush();
      MappedRegion region(std::move(_region));
      _region        = std::move(MappedRegion());
      _isInitialized = false;
    }
  }

private:
  void
  resizeFileAndMapRegion(std::size_t fit_size = 0) {
    std::size_t multiplier = 1;
    std::size_t page_size  = MappedRegion::get_page_size();

    while ((_position + fit_size)
           >= (_chunkSize + multiplier * page_size)) {
      ++multiplier;
    }

    std::lldiv_t divt = std::lldiv(_position, page_size);

    _fileSize     += _chunkSize + multiplier * page_size;
    _position      = divt.rem;
    _regionOffset += divt.quot * page_size;

    FileStream fileStream(_filePath,
                          std::ios_base::in
                          | std::ios_base::out
                          | std::ios_base::binary);

    if (!fileStream) {
      fileStream.close();
      fileStream.open(_filePath,
                      std::fstream::binary
                      | std::fstream::trunc
                      | std::fstream::out);
      fileStream.close();

      fileStream.open(_filePath,
                      std::fstream::binary
                      | std::fstream::in
                      | std::fstream::out);
    }

    std::filebuf* fileBuf = fileStream.rdbuf();
    _fileSize = _fileSize + multiplier * _chunkSize;
    fileBuf->pubseekoff(_fileSize - 1, std::ios_base::beg);
    fileBuf->sputc(0);
    fileStream.close();

    FileMapping _file(_filePath.string().data(),
                      boost::interprocess::read_write);

    MappedRegion region(_file,
                        boost::interprocess::read_write,
                        _regionOffset);
    _region = std::move(region);
  }

  inline std::string
  createRow(std::vector<std::string> const& row) {
    std::string result;

    for (unsigned i = 0; i < row.size(); ++i) {
      if (i != 0) {
        result += ",";
      }
      result += row[i];
    }

    return result;
  }

  inline std::string
  createRow(std::vector < std::pair < std::string, std::string >> const& row) {
    std::string result1;
    std::string result2;

    for (unsigned i = 0; i < row.size(); ++i) {
      if (i != 0) {
        result1 += ",";
        result2 += ",";
      }
      result1 += row[i].first;
      result2 += row[i].second;
    }

    return (result1 + '\n' + result2);
  }

  bool _isInitialized;

  Path         _filePath;
  MappedRegion _region;
  std::size_t  _chunkSize;
  std::size_t  _fileSize;
  std::size_t  _position;
  std::size_t  _regionOffset;

  Path _infoFilePath;

  std::vector<std::string> _row;
  std::vector < std::pair < std::string, std::string >> _info;
};
}
}
