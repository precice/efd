#include "parsefromstring.hpp"

#include <boost/locale.hpp>

using namespace Uni::Helpers;

std::locale const Private::default_locale(
  boost::locale::generator() ("en_US.UTF-8"));

boost::regex const Private::floating_point_number_regex(
  "[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?",
  boost::regbase::ECMAScript);

boost::regex const Private::integer_number_regex(
  "[+-]?\\d+",
  boost::regbase::ECMAScript);
