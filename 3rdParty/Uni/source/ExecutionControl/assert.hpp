#ifndef Uni_ExecutionControl_assert_hpp
#define Uni_ExecutionControl_assert_hpp

#include <Uni/Logging/macros>

namespace Uni {
namespace ExecutionControl {
inline void
AssertionFailed(char const* expr,
                char const* function,
                char const* file,
                int const&  line) {
  Uni_Logging_GLOBAL_SIMPLE_LOG_STREAM(
    Uni_Logging_globalSimpleLogger, 5, file, line)
    << "Expression '"
    << expr << "' failed in function '"
    << function << "'";
  exit(-1);
}
}
}

#ifndef NDEBUG
#  ifndef ASSERT
#    define ASSERT(expr)                                                   \
  if (!expr) {                                                             \
    Uni::ExecutionControl::AssertionFailed(#expr,                          \
                                           __PRETTY_FUNCTION__,            \
                                           Uni_Logging_GLOBAL_LOG_LOCATOR, \
                                           __LINE__)                       \
  }
#    ifndef assert
#      define assert(expr) ASSERT(expr)
#    endif
#  endif
#else
#  ifndef ASSERT
#    define ASSERT(expr)
#    ifndef assert
#      define assert(expr) ASSERT(expr)
#    endif
#  endif
#endif

#ifndef ASSERT_EXIST
#  define ASSERT_EXIST(ptr) ASSERT(ptr != 0)
#  ifndef assertExist
#    define assertExist(expr) ASSERT_EXIST(expr)
#  endif
#endif

#endif
