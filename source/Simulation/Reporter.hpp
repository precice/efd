#pragma once

#include <Eigen/Core>

#include <boost/filesystem.hpp>

namespace Fluid {
namespace Simulation {
class Reporter {
public:
  using Path = boost::filesystem::path;

  Reporter() {}

  Reporter(Reporter const&) = delete;

  virtual ~Reporter() {}

  Reporter&
  operator=(Reporter const&) = delete;

  virtual void
  initialize(Path const& directory_path, std::string file_name_prefix) {
    ((void)directory_path);
    ((void)file_name_prefix);
  }

  template <typename T, int D>
  void
  setAt(std::string const& column_name,
        Eigen::Matrix<T, D, 1> const& value,
        unsigned const& length = 0) {
    this->setAt(column_name, toString(value, length));
  }

  template <typename T>
  void
  setAt(std::string const& column_name,
        T const&           value) {
    this->setAt(column_name, toString(value));
  }

  virtual void
  setAt(std::string const& column_name,
        std::string const& value) {
    ((void)column_name);
    ((void)value);
  }

  template <typename T, int D>
  void
  addAt(std::string const& column,
        Eigen::Matrix<T, D, 1> const& value,
        unsigned const& length = 0) {
    this->addAt(column, toString(value, length));
  }

  template <typename T>
  void
  addAt(std::string const& column, T const& value) {
    this->addAt(column, toString(value));
  }

  virtual void
  addAt(std::string const& column, std::string const& value) {
    ((void)column);
    ((void)value);
  }

  virtual void
  append(std::string const& value) {
    ((void)value);
  }

  virtual void
  recordInfo() {}

  virtual void
  recordIteration() {}

  virtual void
  release() {}

protected:
  template <typename T, int D>
  inline std::string
  toString(Eigen::Matrix<T, D, 1> const& value,
           unsigned const& length = 0) {
    unsigned size = length;

    if (length == 0) {
      size = value.size();
    }

    std::stringstream ss;

    for (unsigned i = 0; i < size; ++i) {
      if (i != 0) {
        ss << " ";
      }
      ss << value(i);
    }

    return ss.str();
  }

  template <typename T>
  inline std::string
  toString(T const& value) {
    std::stringstream ss;
    ss << value;

    return ss.str();
  }
};
}
}
