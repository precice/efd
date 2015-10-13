#ifndef Uni_Logging_macros
#define Uni_Logging_macros

#include <Uni/Logging/LogRecordStreamer>
#include <Uni/Logging/SimpleLogger>

static Uni::Logging::SimpleLogger
  Uni_Logging_globalSimpleLogger("Global");

#define Uni_Logging_GLOBAL_SIMPLE_LOG_STREAM( \
    logger, severity, file_path, line_number) \
  Uni_Logging_LOG_STREAM(                     \
    Uni::Logging::SimpleLogger,               \
    logger, severity, file_path, line_number) \

#define Uni_Logging_GLOBAL_LOG_LOCATOR \
  __FILE__

#define Uni_Logging_GLOBAL_LOG_STREAM(severity) \
  Uni_Logging_GLOBAL_SIMPLE_LOG_STREAM(         \
    Uni_Logging_globalSimpleLogger,             \
    severity, Uni_Logging_GLOBAL_LOG_LOCATOR, __LINE__)

#define Uni_Logging_GLOBAL_LOG(severity, ...)                        \
  Uni_Logging_globalSimpleLogger.log(severity,                       \
                                     Uni_Logging_GLOBAL_LOG_LOCATOR, \
                                     __LINE__,                       \
                                     __VA_ARGS__)

#ifndef DEBUG
#  define DEBUG \
  Uni_Logging_GLOBAL_LOG_STREAM(1)
#endif

#define Uni_Logging_logDebug(...) \
  Uni_Logging_GLOBAL_LOG(1, __VA_ARGS__)

#ifndef logDebug
#  define logDebug(...) \
  Uni_Logging_logDebug(__VA_ARGS__)
#endif

#ifndef INFO
#  define INFO \
  Uni_Logging_GLOBAL_LOG_STREAM(2)
#endif

#define Uni_Logging_logInfo(...) Uni_Logging_GLOBAL_LOG(2, __VA_ARGS__)

#ifndef logInfo
#  define logInfo(...) Uni_Logging_logInfo(__VA_ARGS__)
#endif

#ifndef WARNING
#  define WARNING \
  Uni_Logging_GLOBAL_LOG_STREAM(3)
#endif

#define Uni_Logging_logWarning(...) \
  Uni_Logging_GLOBAL_LOG(3, __VA_ARGS__)

#ifndef logWarning
#  define logWarning(...) \
  Uni_Logging_logWarning(__VA_ARGS__)
#endif

#ifndef ERROR
#  define ERROR                     \
  LOG_STREAM(                       \
    Uni_Logging_globalSimpleLogger, \
    4, Uni_Logging_GLOBAL_LOG_LOCATOR, __LINE__)
#endif

#define Uni_Logging_logError(...) \
  Uni_Logging_GLOBAL_LOG(4, __VA_ARGS__)

#ifndef logError
#  define logError(...) \
  Uni_Logging_logError(__VA_ARGS__)
#endif

#ifndef FATAL
#  define FATAL \
  Uni_Logging_GLOBAL_LOG_STREAM(5)
#endif

#define Uni_Logging_logFatal(...) \
  Uni_Logging_GLOBAL_LOG(5, __VA_ARGS__)

#ifndef logFatal
#  define logFatal(...) \
  Uni_Logging_logFatal(__VA_ARGS__)
#endif

#endif
