#ifndef Uni_Logging_SimpleLogger
#define Uni_Logging_SimpleLogger

#include "format.hpp"

#include <iostream>

namespace Uni {
namespace Logging {
class SimpleLogger {
public:
  SimpleLogger(char const* name) : _name(name) {}

  template <typename... TArgs>
  void
  log(int const&         severity,
      char const*        notifiedFilePath,
      int const&         notifiedLineNumber,
      std::string const& message,
      TArgs...           args) {
    std::string severity_message = "Unknown";

    switch (severity) {
    case (1):
      severity_message = "Debug ";
      break;

    case (2):
      severity_message = "Info  ";
      break;

    case (3):
      severity_message = "Warning";
      break;

    case (4):
      severity_message = "Error  ";
      break;

    case (5):
      severity_message = "Fatal  ";
      break;
    }

    auto headerStr = format("[{1} | {2}]: ({3}, {4})",
                            _name,
                            severity_message,
                            notifiedFilePath,
                            notifiedLineNumber);

    auto messageStr = format(message, args...);

    std::cout << headerStr.c_str() << std::endl
              << messageStr.c_str() << std::endl << std::flush;
  }

private:
  std::string _name;
};
}
}

#endif
