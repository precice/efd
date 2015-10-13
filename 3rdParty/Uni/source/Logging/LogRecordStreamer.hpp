#ifndef Uni_LogRecordStreamer_LogRecordStreamer
#define Uni_LogRecordStreamer_LogRecordStreamer

#include <sstream>

namespace Uni {
namespace Logging {
template <typename Logger>
class LogRecordStreamer {
public:
  LogRecordStreamer(Logger&     logger,
                    int const&  severity,
                    char const* file_path,
                    int const&  line_number) : _is_first_time(true),
                                               _logger(logger),
                                               _severity(severity),
                                               _file_path(file_path),
                                               _line_number(line_number) {}

  ~LogRecordStreamer() {
    flush();
  }

  void
  flush() {
    _logger.log(_severity,
                _file_path,
                _line_number,
                _stream.str());
    _stream.flush();
  }

  operator bool() {
    if (_is_first_time) {
      _is_first_time = false;

      return true;
    }

    return false;
  }

  template <typename T>
  LogRecordStreamer&
  operator<<(T const& object) {
    _stream << object;

    return *this;
  }

private:
  bool               _is_first_time;
  Logger&            _logger;
  int const          _severity;
  char const*        _file_path;
  int const          _line_number;
  std::ostringstream _stream;
};
}
}

#define Uni_Logging_LOG_STREAM(                       \
    Logger, logger, severity, file_path, line_number) \
  for (Uni::Logging::LogRecordStreamer<Logger>        \
       Uni_Logging_stream_(logger,                    \
                           severity,                  \
                           file_path,                 \
                           line_number);              \
       Uni_Logging_stream_;)                          \
    Uni_Logging_stream_

#endif
// vim:ft=cpp:fenc=utf-8:ff=unix:ts=2:sw=2:tw=80:et:
