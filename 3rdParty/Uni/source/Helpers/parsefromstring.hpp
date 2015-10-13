#pragma once

#include <Eigen/Core>

#include <boost/regex.hpp>

#include <string>
#include <sstream>

namespace Uni {
namespace Helpers {
namespace Private {
extern std::locale const  default_locale;
extern boost::regex const floating_point_number_regex;
extern boost::regex const integer_number_regex;
}

template <typename T>
bool
parse_integer_number(std::string const& str,
                    T&                 result) {
  boost::smatch match;

  if (boost::regex_search(str,
                          match,
                          Private::integer_number_regex)) {
    std::stringstream ss;
    ss.imbue(Private::default_locale);
    ss << match.str();
    ss >> result;

    return true;
  }

  return false;
}

template <typename T>
bool
parse_floating_point_number(std::string const& str,
                            T&                 result) {
  boost::smatch match;

  if (boost::regex_search(str,
                          match,
                          Private::floating_point_number_regex)) {
    std::stringstream ss;
    ss.imbue(Private::default_locale);
    ss << match.str();
    ss >> result;

    return true;
  }

  return false;
}

template <typename T, int D>
int
parse_integer_vector(
  std::string const& str,
  Eigen::Matrix<T, D, 1>&                 result) {
  auto it = boost::sregex_iterator(str.begin(),
                                   str.end(),
                                   Private::integer_number_regex);

  auto it_end = boost::sregex_iterator();

  int i = 0;

  for (; it != it_end; ++it) {
    if (i >= D) {
      return i;
    }
    auto const& match = *it;

    std::stringstream ss;
    ss.imbue(Private::default_locale);
    ss << match.str();
    ss >> result(i);
    ++i;
  }

  return i;
}

template <typename T, int D>
int
parse_floating_point_vector(
  std::string const& str,
  Eigen::Matrix<T, D, 1>&                 result) {
  auto it = boost::sregex_iterator(str.begin(),
                                   str.end(),
                                   Private::floating_point_number_regex);

  auto it_end = boost::sregex_iterator();

  int i = 0;

  for (; it != it_end; ++it) {
    if (i >= D) {
      return i;
    }
    auto const& match = *it;

    std::stringstream ss;
    ss.imbue(Private::default_locale);
    ss << match.str();
    ss >> result(i);
    ++i;
  }

  return i;
}
}
}
