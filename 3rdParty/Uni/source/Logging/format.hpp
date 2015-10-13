#pragma once

#include <boost/locale.hpp>

namespace Uni {
namespace Logging {
namespace Private {
class Formatter {
private:
  typedef boost::locale::generator LocaleGenerator;
  typedef boost::locale::format    Format;

public:
  template <typename... TArgs>
  static inline std::string
  format(std::string const& message,
         TArgs&& ...        args) {
    return get()._format(message, std::forward<TArgs>(args) ...);
  }

private:
  template <typename... TArgs>
  inline std::string
  _format(std::string const& message,
          TArgs&& ...        args) {
    Format messageFormat(message);
    _format(messageFormat, std::forward<TArgs>(args) ...);

    return messageFormat.str(_locale);
  }

  Formatter() : _localeGenerator(),
                _locale(_localeGenerator("en_US.UTF-8")) {}

  Formatter(Formatter const& other) = delete;

  ~Formatter() {}

  Formatter&
  operator=(Formatter const& other) = delete;

  template <typename TValue, typename... TArgs>
  void
  _format(Format& format, TValue const& arg, TArgs&& ... args) {
    format % arg;
    _format(format, std::forward<TArgs>(args) ...);
  }

  void
  _format(Format& format) {
    ((void)format);
  }

  static Formatter&
  get() {
    static Formatter _instance;

    return _instance;
  }

  LocaleGenerator _localeGenerator;
  std::locale     _locale;
};
}

template <typename... TArgs>
std::string
format(std::string const& message,
       TArgs&& ...        args) {
  return Private::Formatter::format(message, std::forward<TArgs>(args) ...);
}
}
}
