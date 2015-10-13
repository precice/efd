#ifndef Uni_ExecutionControl_exception_hpp
#define Uni_ExecutionControl_exception_hpp

#include <Uni/Logging/macros>

#include <exception>

namespace Uni {
namespace ExecutionControl {
inline void
throwException() {
  throw std::exception();
}
}
}

#define Uni_ExecutionControl_EXCEPTION_STREAM(                         \
    Logger, logger, severity, file_path, line_number)                  \
  for (Uni::Logging::LogRecordStreamer<Logger>                         \
       Uni_Auxiliary_stream_(logger,                                   \
                             severity,                                 \
                             file_path,                                \
                             line_number);                             \
       Uni_Auxiliary_stream_; Uni::ExecutionControl::throwException()) \
    Uni_Auxiliary_stream_

#ifndef EXCEPTION
#  define EXCEPTION                      \
  Uni_ExecutionControl_EXCEPTION_STREAM( \
    Uni::Logging::SimpleLogger,          \
    _Uni_Logging_globalSimpleLogger,     \
    5, Uni_Logging_GLOBAL_LOG_LOCATOR, __LINE__)

#endif

#ifndef throwException
#  define throwException(...)             \
  Uni_Logging_GLOBAL_LOG(1, __VA_ARGS__); \
  Uni::ExecutionControl::throwException()
#endif

#endif
