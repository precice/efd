#ifndef Uni_Platform_compiler_hpp
#define Uni_Platform_compiler_hpp

#if                                                                            \
  defined (__MINGW32__)                                                     || \
  defined (__MINGW64__)
#  define Uni_Platform_COMPILER "MinGW"
#  define Uni_Platform_COMPILER_MINGW
#elif \
  defined (__GNUC__)
#  define Uni_Platform_COMPILER "GNU"
#  define Uni_Platform_COMPILER_GNU
#  define Uni_Platform_COMPILER_GNU_VERSION_MAJOR __GNUC__
#  define Uni_Platform_COMPILER_GNU_VERSION_MINOR __GNUC_MINOR__
#  define Uni_Platform_COMPILER_GNU_VERSION_PATCH __GNUC_PATCHLEVEL__
#elif \
  defined (__clang__)
#  define Uni_Platform_COMPILER "Clang"
#  define Uni_Platform_COMPILER_CLANG
#elif \
  defined (_MSC_VER)
#  define Uni_Platform_COMPILER "Microsoft Visual C++"
#  define Uni_Platform_COMPILER_MICROSOFT
#elif \
  defined (__BORLANDC__)
#  define Uni_Platform_COMPILER "Borland C++ Builder"
#  define Uni_Platform_COMPILER_BORLAND
#elif \
  defined (__CODEGEARC__)
#  define Uni_Platform_COMPILER "CodeGear C++ Builder"
#  define Uni_Platform_COMPILER_CODEGEAR
#elif                                                                          \
  defined (__INTEL_COMPILER)                                                || \
  defined (__ICL)
#  define Uni_Platform_COMPILER "Intel C++"
#  define Uni_Platform_COMPILER_INTEL
#elif                                                                          \
  defined (__xlC__)                                                         || \
  defined (__IBMCPP__)
#  define Uni_Platform_COMPILER "IBM XL C++"
#  define Uni_Platform_COMPILER_IBM
#elif \
  defined (__HP_aCC)
#  define Uni_Platform_COMPILER "HP aC++"
#  define Uni_Platform_COMPILER_HP
#elif \
  defined (__WATCOMC__)
#  define Uni_Platform_COMPILER "Watcom C++"
#  define Uni_Platform_COMPILER_WATCOM
#endif

#ifndef Uni_Platform_COMPILER
#  error "Current compiler is not supported."
#endif

#endif
