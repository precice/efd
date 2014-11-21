#ifndef FsiSimulation_EntryPoint_Private_convertUtfPathToAnsi_hpp
#define FsiSimulation_EntryPoint_Private_convertUtfPathToAnsi_hpp

#include <boost/filesystem.hpp>
#include <boost/locale.hpp>

namespace FsiSimulation {
namespace EntryPoint {
namespace Private {
inline std::string
convertUtfPathToAnsi(std::string const& input,
                     std::locale const& utfLocale,
                     std::locale const& ansiLocale) {
  auto fromEncoding = std::use_facet<boost::locale::info>(
    utfLocale).encoding();
  auto toEncoding = std::use_facet<boost::locale::info>(
    ansiLocale).encoding();

  return boost::locale::conv::between(input,
                                      toEncoding, fromEncoding);
}

inline std::string
convertUtfPathToAnsi(std::string const& input) {
  using LocaleGenerator = boost::locale::generator;

  auto globalLocale = LocaleGenerator().generate("");

  LocaleGenerator ansiLocaleGenerator;
  ansiLocaleGenerator.use_ansi_encoding(true);
  auto ansiLocale = ansiLocaleGenerator.generate("");

  return convertUtfPathToAnsi(input,
                              globalLocale,
                              ansiLocale);
}
}
}
}

namespace boost {
namespace filesystem {
inline boost::filesystem::path
make_relative(boost::filesystem::path a_To,
              boost::filesystem::path a_From =
                boost::filesystem::   current_path()) {
  a_To   = boost::filesystem::canonical(a_To);
  a_From = boost::filesystem::canonical(a_From);
  boost::filesystem::path                 ret;
  boost::filesystem::path::const_iterator itrFrom(a_From.begin());
  boost::filesystem::path::const_iterator itrTo(a_To.begin());

  // Find common base
  for (boost::filesystem::path::const_iterator toEnd(a_To.end()),
       fromEnd(a_From.end());
       itrFrom != fromEnd && itrTo != toEnd && *itrFrom == *itrTo;
       ++itrFrom, ++itrTo) {}

  // Navigate backwards in directory to reach previously found base
  for (boost::filesystem::path::const_iterator fromEnd(a_From.end());
       itrFrom != fromEnd;
       ++itrFrom) {
    if ((*itrFrom) != ".") {
      ret /= "..";
    }
  }

  // Now navigate down the directory branch
  for (; itrTo != a_To.end(); ++itrTo) {
    ret /= *itrTo;
  }

  return ret;
}
}
}

#endif
