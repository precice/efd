#ifndef Uni_Platform_operatingsystem_hpp
#define Uni_Platform_operatingsystem_hpp

#if                                                                            \
  defined (__CYGWIN__)                                                      || \
  defined (__CYGWIN32__)
#  define Uni_Platform_OS "Cygwin"
#  define Uni_Platform_OS_CYGWIN
#  define Uni_Platform_OS_UNIX
#  define Uni_Platform_OS_WINDOWS
#elif                                                                          \
  defined (_WIN16)                                                          || \
  defined (_WIN32)                                                          || \
  defined (_WIN64)                                                          || \
  defined (__WIN32__)                                                       || \
  defined (__TOS_WIN__)                                                     || \
  defined (__WINDOWS__)
#  define Uni_Platform_OS "Windows"
#  define Uni_Platform_OS_WINDOWS
#elif                                                                          \
  defined (macintosh)                                                       || \
  defined (Macintosh)                                                       || \
  defined (__TOS_MACOS__)                                                   || \
  (defined (__APPLE__) && defined (__MACH__))
#  define Uni_Platform_OS "Mac"
#  define Uni_Platform_OS_MAC
#  define Uni_Platform_OS_UNIX
#elif                                                                          \
  defined (linux)                                                           || \
  defined (__linux)                                                         || \
  defined (__linux__)                                                       || \
  defined (__TOS_LINUX__)
#  define Uni_Platform_OS "Linux"
#  define Uni_Platform_OS_LINUX
#  define Uni_Platform_OS_UNIX
#elif                                                                          \
  defined (__FreeBSD__)                                                     || \
  defined (__OpenBSD__)                                                     || \
  defined (__NetBSD__)                                                      || \
  defined (__bsdi__)                                                        || \
  defined (__DragonFly__)
#  define Uni_Platform_OS "BSD"
#  define Uni_Platform_OS_BSD
#  define Uni_Platform_OS_UNIX
#elif                                                                          \
  defined (sun)                                                             || \
  defined (__sun)
#  define Uni_Platform_OS "Solaris"
#  define Uni_Platform_OS_SOLARIS
#  define Uni_Platform_OS_UNIX
#elif                                                                          \
  defined (_AIX)                                                            || \
  defined (__TOS_AIX__)
#  define Uni_Platform_OS "AIX"
#  define Uni_Platform_OS_AIX
#  define Uni_Platform_OS_UNIX
#elif                                                                          \
  defined (hpux)                                                            || \
  defined (_hpux)                                                           || \
  defined (__hpux)
#  define Uni_Platform_OS "HPUX"
#  define Uni_Platform_OS_HPUX
#  define Uni_Platform_OS_UNIX
#elif \
  defined (__QNX__)
#  define Uni_Platform_OS "QNX"
#  define Uni_Platform_OS_QNX
#  define Uni_Platform_OS_UNIX
#elif                                                                          \
  defined (unix)                                                            || \
  defined (__unix)                                                          || \
  defined (__unix__)
#  define Uni_Platform_OS "Unix"
#  define Uni_Platform_OS_UNIX
#endif

#ifndef Uni_Platform_OS
#  error "Current platform is not supported."
#endif

#endif
