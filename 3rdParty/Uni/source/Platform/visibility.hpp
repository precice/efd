#ifndef Uni_Platform_visibility_hpp
#define Uni_Platform_visibility_hpp

#include <Uni/Platform/compiler>
#include <Uni/Platform/operatingsystem>

#ifdef Uni_Platform_OS_WINDOWS
#  define Uni_Platform_EXPORT __declspec(dllexport)
#  define Uni_Platform_IMPORT __declspec(dllimport)
#  define Uni_Platform_LOCAL
#elif                                                                     \
  Uni_Platform_COMPILER_GNU_VERSION_MAJOR >= 4                         || \
  defined (Uni_Platform_COMPILER_CLANG)
#  define Uni_Platform_EXPORT __attribute__((visibility("default")))
#  define Uni_Platform_IMPORT __attribute__((visibility("default")))
#  define Uni_Platform_LOCAL  __attribute__((visibility("hidden")))
#else
#  define Uni_Platform_EXPORT
#  define Uni_Platform_IMPORT
#  define Uni_Platform_LOCAL
#endif

// #define Uni_Platform_EXPORT __declspec(dllexport)
// #define Uni_Platform_IMPORT __declspec(dllimport)
// #define Uni_Platform_LOCAL

#ifdef __cplusplus
#  define Uni_Platform_DEMANGLED extern "C"
#else
#  define Uni_Platform_DEMANGLED
#endif

#endif
