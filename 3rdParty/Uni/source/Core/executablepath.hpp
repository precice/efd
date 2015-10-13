#pragma once

#include <Uni/ExecutionControl/exception>
#include <Uni/Platform/operatingsystem>
#include <boost/filesystem.hpp>

#if defined (Uni_Platform_OS_WINDOWS)
#  include <windows.h>
#elif defined (Uni_Platform_OS_UNIX)

#else
#  error "Unknown platform"
#endif

//

namespace Uni {
#if defined (Uni_Platform_OS_WINDOWS)
inline boost::filesystem::path
getExecutablePath() {
  TCHAR buffer[MAX_PATH];
  auto  size =  GetModuleFileNameW(NULL, buffer, MAX_PATH);

  if (size == MAX_PATH || size == 0) {
    throwException("Failed to locate application");
  }
  namespace fs = boost::filesystem;
  fs::path result(buffer);

  return result;
}
#elif defined (Uni_Platform_OS_UNIX)
inline boost::filesystem::path
getExecutablePath() {
  char buffer[PATH_MAX];
  auto size = readlink("/proc/self/exe", buffer, PATH_MAX);

  if (size == -1) {
    throwException("Failed to locate application");
  }

  buffer[size] = '\0';

  namespace fs = boost::filesystem;
  fs::path result(buffer);

  return result;
}

#else
#  error "Unknown platform"
#endif
}
