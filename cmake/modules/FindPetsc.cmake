# - Try to find Petsc
# Once done this will define
#
#  PETSC_FOUND        - system has Petsc
#  PETSC_INCLUDES     - the Petsc include directories
#  PETSC_LIBRARIES    - Link these to use Petsc
#  PETSC_COMPILER     - Compiler used by Petsc, helpful to find a compatible MPI
#  PETSC_DEFINITIONS  - Compiler switches for using Petsc
#  PETSC_MPIEXEC      - Executable for running MPI programs
#  PETSC_VERSION      - Version string (MAJOR.MINOR.SUBMINOR)
#
#  Usage:
#  find_package(Petsc COMPONENTS CXX)  - required if build --with-clanguage=C++ --with-c-support=0
#  find_package(Petsc COMPONENTS C)    - standard behavior of checking build using a C compiler
#  find_package(Petsc)                 - same as above
#
# Setting these changes the behavior of the search
#  PETSC_DIR - directory in which Petsc resides
#  PETSC_ARCH - build architecture
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file. #


# FindPackageMultipass {{{
# ============================================================================
# PackageMultipass - this module defines two macros
#
# FIND_PACKAGE_MULTIPASS (Name CURRENT
#  STATES VAR0 VAR1 ...
#  DEPENDENTS DEP0 DEP1 ...)
#
#  This function creates a cache entry <UPPERCASED-Name>_CURRENT which
#  the user can set to "NO" to trigger a reconfiguration of the package.
#  The first time this function is called, the values of
#  <UPPERCASED-Name>_VAR0, ... are saved.  If <UPPERCASED-Name>_CURRENT
#  is false or if any STATE has changed since the last time
#  FIND_PACKAGE_MULTIPASS() was called, then CURRENT will be set to "NO",
#  otherwise CURRENT will be "YES".  IF not CURRENT, then
#  <UPPERCASED-Name>_DEP0, ... will be FORCED to NOTFOUND.
#  Example:
#    find_path (FOO_DIR include/foo.h)
#    FIND_PACKAGE_MULTIPASS (Foo foo_current
#      STATES DIR
#      DEPENDENTS INCLUDES LIBRARIES)
#    if (NOT foo_current)
#      # Make temporary files, run programs, etc, to determine FOO_INCLUDES and FOO_LIBRARIES
#    endif (NOT foo_current)
#
# MULTIPASS_SOURCE_RUNS (Name INCLUDES LIBRARIES SOURCE RUNS LANGUAGE)
#  Always runs the given test, use this when you need to re-run tests
#  because parent variables have made old cache entries stale. The LANGUAGE
#  variable is either C or CXX indicating which compiler the test should
#  use.
# MULTIPASS_C_SOURCE_RUNS (Name INCLUDES LIBRARIES SOURCE RUNS)
#  DEPRECATED! This is only included for backwards compatability. Use
#  the more general MULTIPASS_SOURCE_RUNS instead.
#  Always runs the given test, use this when you need to re-run tests
#  because parent variables have made old cache entries stale.

macro (FIND_PACKAGE_MULTIPASS _name _current)
  string (TOUPPER ${_name} _NAME)
  set (_args ${ARGV})
  list (REMOVE_AT _args 0 1)

  set (_states_current "YES")
  list (GET _args 0 _cmd)
  if (_cmd STREQUAL "STATES")
    list (REMOVE_AT _args 0)
    list (GET _args 0 _state)
    while (_state AND NOT _state STREQUAL "DEPENDENTS")
      # The name of the stored value for the given state
      set (_stored_var PACKAGE_MULTIPASS_${_NAME}_${_state})
      if (NOT "${${_stored_var}}" STREQUAL "${${_NAME}_${_state}}")
        set (_states_current "NO")
      endif (NOT "${${_stored_var}}" STREQUAL "${${_NAME}_${_state}}")
      set (${_stored_var} "${${_NAME}_${_state}}" CACHE INTERNAL "Stored state for ${_name}." FORCE)
      list (REMOVE_AT _args 0)
      list (GET _args 0 _state)
    endwhile (_state AND NOT _state STREQUAL "DEPENDENTS")
  endif (_cmd STREQUAL "STATES")

  set (_stored ${_NAME}_CURRENT)
  if (NOT ${_stored})
    set (${_stored} "YES" CACHE BOOL "Is the configuration for ${_name} current?  Set to \"NO\" to reconfigure." FORCE)
    set (_states_current "NO")
  endif (NOT ${_stored})

  set (${_current} ${_states_current})
  if (NOT ${_current} AND PACKAGE_MULTIPASS_${_name}_CALLED)
    message (STATUS "Clearing ${_name} dependent variables")
    # Clear all the dependent variables so that the module can reset them
    list (GET _args 0 _cmd)
    if (_cmd STREQUAL "DEPENDENTS")
      list (REMOVE_AT _args 0)
      foreach (dep ${_args})
        set (${_NAME}_${dep} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
      endforeach (dep)
    endif (_cmd STREQUAL "DEPENDENTS")
    set (${_NAME}_FOUND "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
  endif ()
  set (PACKAGE_MULTIPASS_${name}_CALLED YES CACHE INTERNAL "Private" FORCE)
endmacro (FIND_PACKAGE_MULTIPASS)


macro (MULTIPASS_SOURCE_RUNS includes libraries source runs language)
  include (Check${language}SourceRuns)
  # This is a ridiculous hack.  CHECK_${language}_SOURCE_* thinks that if the
  # *name* of the return variable doesn't change, then the test does
  # not need to be re-run.  We keep an internal count which we
  # increment to guarantee that every test name is unique.  If we've
  # gotten here, then the configuration has changed enough that the
  # test *needs* to be rerun.
  if (NOT MULTIPASS_TEST_COUNT)
    set (MULTIPASS_TEST_COUNT 00)
  endif (NOT MULTIPASS_TEST_COUNT)
  math (EXPR _tmp "${MULTIPASS_TEST_COUNT} + 1") # Why can't I add to a cache variable?
  set (MULTIPASS_TEST_COUNT ${_tmp} CACHE INTERNAL "Unique test ID")
  set (testname MULTIPASS_TEST_${MULTIPASS_TEST_COUNT}_${runs})
  set (CMAKE_REQUIRED_INCLUDES ${includes})
  set (CMAKE_REQUIRED_LIBRARIES ${libraries})
  if(${language} STREQUAL "C")
    check_c_source_runs ("${source}" ${testname})
  elseif(${language} STREQUAL "CXX")
    check_cxx_source_runs ("${source}" ${testname})
  endif()
  set (${runs} "${${testname}}")
endmacro (MULTIPASS_SOURCE_RUNS)

macro (MULTIPASS_C_SOURCE_RUNS includes libraries source runs)
  multipass_source_runs("${includes}" "${libraries}" "${source}" ${runs} "C")
endmacro (MULTIPASS_C_SOURCE_RUNS)
# ============================================================================
# }}} FindPackageMultipass

# CorrectWindowsPaths {{{
# ============================================================================
# CorrectWindowsPaths - this module defines one macro
#
# CONVERT_CYGWIN_PATH( PATH )
# This uses the command cygpath (provided by cygwin) to convert
# unix-style paths into paths useable by cmake on windows
macro (CONVERT_CYGWIN_PATH _path)
  if (WIN32)
    EXECUTE_PROCESS(COMMAND cygpath.exe -m ${${_path}}
    OUTPUT_VARIABLE ${_path})
    string (STRIP ${${_path}} ${_path})
  endif (WIN32)
endmacro (CONVERT_CYGWIN_PATH)
# ============================================================================
# }}} CorrectWindowsPaths

# ResolveCompilerPaths {{{
# ============================================================================
# ResolveCompilerPaths - this module defines two macros
#
# RESOLVE_LIBRARIES (XXX_LIBRARIES LINK_LINE)
#  This macro is intended to be used by FindXXX.cmake modules.
#  It parses a compiler link line and resolves all libraries
#  (-lfoo) using the library path contexts (-L/path) in scope.
#  The result in XXX_LIBRARIES is the list of fully resolved libs.
#  Example:
#
#    RESOLVE_LIBRARIES (FOO_LIBRARIES "-L/A -la -L/B -lb -lc -ld")
#
#  will be resolved to
#
#    FOO_LIBRARIES:STRING="/A/liba.so;/B/libb.so;/A/libc.so;/usr/lib/libd.so"
#
#  if the filesystem looks like
#
#    /A:       liba.so         libc.so
#    /B:       liba.so libb.so
#    /usr/lib: liba.so libb.so libc.so libd.so
#
#  and /usr/lib is a system directory.
#
#  Note: If RESOLVE_LIBRARIES() resolves a link line differently from
#  the native linker, there is a bug in this macro (please report it).
#
# RESOLVE_INCLUDES (XXX_INCLUDES INCLUDE_LINE)
#  This macro is intended to be used by FindXXX.cmake modules.
#  It parses a compile line and resolves all includes
#  (-I/path/to/include) to a list of directories.  Other flags are ignored.
#  Example:
#
#    RESOLVE_INCLUDES (FOO_INCLUDES "-I/A -DBAR='\"irrelevant -I/string here\"' -I/B")
#
#  will be resolved to
#
#    FOO_INCLUDES:STRING="/A;/B"
#
#  assuming both directories exist.
#  Note: as currently implemented, the -I/string will be picked up mistakenly (cry, cry)
macro (RESOLVE_LIBRARIES LIBS LINK_LINE)
  string (REGEX MATCHALL "((-L|-l|-Wl)([^\" ]+|\"[^\"]+\")|[^\" ]+\\.(a|so|dll|lib))" _all_tokens "${LINK_LINE}")
  set (_libs_found)
  set (_directory_list)
  foreach (token ${_all_tokens})
    if (token MATCHES "-L([^\" ]+|\"[^\"]+\")")
      # If it's a library path, add it to the list
      string (REGEX REPLACE "^-L" "" token ${token})
      string (REGEX REPLACE "//" "/" token ${token})
      convert_cygwin_path(token)
      list (APPEND _directory_list ${token})
    elseif (token MATCHES "^(-l([^\" ]+|\"[^\"]+\")|[^\" ]+\\.(a|so|dll|lib))")
      # It's a library, resolve the path by looking in the list and then (by default) in system directories
      if (WIN32) #windows expects "libfoo", linux expects "foo"
        string (REGEX REPLACE "^-l" "lib" token ${token})
      else (WIN32)
        string (REGEX REPLACE "^-l" "" token ${token})
      endif (WIN32)
      set (_root)
      if (token MATCHES "^/")	# We have an absolute path
        #separate into a path and a library name:
        string (REGEX MATCH "[^/]*\\.(a|so|dll|lib)$" libname ${token})
        string (REGEX MATCH ".*[^${libname}$]" libpath ${token})
        convert_cygwin_path(libpath)
        set (_directory_list ${_directory_list} ${libpath})
        set (token ${libname})
      endif (token MATCHES "^/")
      set (_lib "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
      find_library (_lib ${token} HINTS ${_directory_list} ${_root})
      if (_lib)
	string (REPLACE "//" "/" _lib ${_lib})
        list (APPEND _libs_found ${_lib})
      else (_lib)
        message (STATUS "Unable to find library ${token}")
      endif (_lib)
    endif (token MATCHES "-L([^\" ]+|\"[^\"]+\")")
  endforeach (token)
  set (_lib "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)
  # only the LAST occurence of each library is required since there should be no circular dependencies
  if (_libs_found)
    list (REVERSE _libs_found)
    list (REMOVE_DUPLICATES _libs_found)
    list (REVERSE _libs_found)
  endif (_libs_found)
  set (${LIBS} "${_libs_found}")
endmacro (RESOLVE_LIBRARIES)

macro (RESOLVE_INCLUDES INCS COMPILE_LINE)
  string (REGEX MATCHALL "-I([^\" ]+|\"[^\"]+\")" _all_tokens "${COMPILE_LINE}")
  set (_incs_found)
  foreach (token ${_all_tokens})
    string (REGEX REPLACE "^-I" "" token ${token})
    string (REGEX REPLACE "//" "/" token ${token})
    convert_cygwin_path(token)
    if (EXISTS ${token})
      list (APPEND _incs_found ${token})
    else (EXISTS ${token})
      message (STATUS "Include directory ${token} does not exist")
    endif (EXISTS ${token})
  endforeach (token)
  list (REMOVE_DUPLICATES _incs_found)
  set (${INCS} "${_incs_found}")
endmacro (RESOLVE_INCLUDES)
# ============================================================================
# }}} ResolveCompilerPaths

set(PETSC_VALID_COMPONENTS
    C
    CXX
    )

if(NOT Petsc_FIND_COMPONENTS)
  set(PETSC_LANGUAGE_BINDINGS "C")
else()
  # Right now, this is designed for compatability with the --with-clanguage option, so
  # only allow one item in the components list.
  list(LENGTH ${Petsc_FIND_COMPONENTS} components_length)
  if(${components_length} GREATER 1)
    message(FATAL_ERROR "Only one component for Petsc is allowed to be specified")
  endif()
  # This is a stub for allowing multiple components should that time ever come. Perhaps
  # to also test Fortran bindings?
  foreach(component ${Petsc_FIND_COMPONENTS})
    list(FIND PETSC_VALID_COMPONENTS ${component} component_location)
    if(${component_location} EQUAL -1)
      message(FATAL_ERROR "\"${component}\" is not a valid Petsc component.")
    else()
      list(APPEND PETSC_LANGUAGE_BINDINGS ${component})
    endif()
  endforeach()
endif()

function (petsc_get_version)
  if (EXISTS "${PETSC_DIR}/include/petscversion.h")
    file (STRINGS "${PETSC_DIR}/include/petscversion.h" vstrings REGEX "#define PETSC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
    foreach (line ${vstrings})
      string (REGEX REPLACE " +" ";" fields ${line}) # break line into three fields (the first is always "#define")
      list (GET fields 1 var)
      list (GET fields 2 val)
      set (${var} ${val} PARENT_SCOPE)
      set (${var} ${val})         # Also in local scope so we have access below
    endforeach ()
    if (PETSC_VERSION_RELEASE)
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}p${PETSC_VERSION_PATCH}" PARENT_SCOPE)
    else ()
      # make dev version compare higher than any patch level of a released version
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.99" PARENT_SCOPE)
    endif ()
  else ()
    message (SEND_ERROR "PETSC_DIR can not be used, ${PETSC_DIR}/include/petscversion.h does not exist")
  endif ()
endfunction ()

find_path (PETSC_DIR include/petsc.h
  HINTS ENV PETSC_DIR
  PATHS
  # Debian paths
  /usr/lib/petscdir/3.5.1 /usr/lib/petscdir/3.5
  /usr/lib/petscdir/3.4.2 /usr/lib/petscdir/3.4
  /usr/lib/petscdir/3.3 /usr/lib/petscdir/3.2 /usr/lib/petscdir/3.1
  /usr/lib/petscdir/3.0.0 /usr/lib/petscdir/2.3.3 /usr/lib/petscdir/2.3.2
  # MacPorts path
  /opt/local/lib/petsc
  $ENV{HOME}/petsc
  DOC "Petsc Directory")

find_program (MAKE_EXECUTABLE NAMES make gmake)

if (PETSC_DIR AND NOT PETSC_ARCH)
  set (_petsc_arches
    $ENV{PETSC_ARCH}                   # If set, use environment variable first
    linux-gnu-c-debug linux-gnu-c-opt  # Debian defaults
    x86_64-unknown-linux-gnu i386-unknown-linux-gnu)
  set (petscconf "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
  foreach (arch ${_petsc_arches})
    if (NOT PETSC_ARCH)
      find_path (petscconf petscconf.h
        HINTS ${PETSC_DIR}
        PATH_SUFFIXES ${arch}/include bmake/${arch}
        NO_DEFAULT_PATH)
      if (petscconf)
        set (PETSC_ARCH "${arch}" CACHE STRING "Petsc build architecture")
      endif (petscconf)
    endif (NOT PETSC_ARCH)
  endforeach (arch)
  set (petscconf "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)
endif (PETSC_DIR AND NOT PETSC_ARCH)

set (petsc_slaves LIBRARIES_SYS LIBRARIES_VEC LIBRARIES_MAT LIBRARIES_DM LIBRARIES_KSP LIBRARIES_SNES LIBRARIES_TS
  INCLUDE_DIR INCLUDE_CONF)
find_package_multipass (Petsc petsc_config_current
  STATES DIR ARCH
  DEPENDENTS INCLUDES LIBRARIES COMPILER MPIEXEC ${petsc_slaves})

# Determine whether the Petsc layout is old-style (through 2.3.3) or
# new-style (>= 3.0.0)
if (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h")   # > 2.3.3
  set (petsc_conf_rules "${PETSC_DIR}/conf/rules")
  set (petsc_conf_variables "${PETSC_DIR}/conf/variables")
elseif (EXISTS "${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf.h") # <= 2.3.3
  set (petsc_conf_rules "${PETSC_DIR}/bmake/common/rules")
  set (petsc_conf_variables "${PETSC_DIR}/bmake/common/variables")
elseif (PETSC_DIR)
  message (SEND_ERROR "The pair PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} do not specify a valid Petsc installation")
endif ()

if (petsc_conf_rules AND petsc_conf_variables AND NOT petsc_config_current)
  petsc_get_version()

  # Put variables into environment since they are needed to get
  # configuration (petscvariables) in the Petsc makefile
  set (ENV{PETSC_DIR} "${PETSC_DIR}")
  set (ENV{PETSC_ARCH} "${PETSC_ARCH}")

  # A temporary makefile to probe the Petsc configuration
  set (petsc_config_makefile "${PROJECT_BINARY_DIR}/Makefile.petsc")
  file (WRITE "${petsc_config_makefile}"
"## This file was autogenerated by FindPetsc.cmake
# PETSC_DIR  = ${PETSC_DIR}
# PETSC_ARCH = ${PETSC_ARCH}
include ${petsc_conf_rules}
include ${petsc_conf_variables}
show :
\t-@echo -n \${\${VARIABLE}}
")

  macro (PETSC_GET_VARIABLE name var)
    set (${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
    execute_process (COMMAND ${MAKE_EXECUTABLE} --no-print-directory -f ${petsc_config_makefile} show VARIABLE=${name}
      OUTPUT_VARIABLE ${var}
      RESULT_VARIABLE petsc_return)
  endmacro (PETSC_GET_VARIABLE)
  petsc_get_variable (PETSC_LIB_DIR            petsc_lib_dir)
  petsc_get_variable (PETSC_EXTERNAL_LIB_BASIC petsc_libs_external)
  petsc_get_variable (PETSC_CCPPFLAGS          petsc_cpp_line)
  petsc_get_variable (PETSC_INCLUDE            petsc_include)
  petsc_get_variable (PCC                      petsc_cc)
  petsc_get_variable (PCC_FLAGS                petsc_cc_flags)
  petsc_get_variable (MPIEXEC                  petsc_mpiexec)
  # We are done with the temporary Makefile, calling PETSC_GET_VARIABLE after this point is invalid!
  file (REMOVE ${petsc_config_makefile})

  # Extract include paths and libraries from compile command line
  resolve_includes (petsc_includes_all "${petsc_cpp_line}")

  #on windows we need to make sure we're linking against the right
  #runtime library
  if (WIN32)
    if (petsc_cc_flags MATCHES "-MT")
      set(using_md False)
      foreach(flag_var
          CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
          CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        if(${flag_var} MATCHES "/MD")
          set(using_md True)
        endif(${flag_var} MATCHES "/MD")
      endforeach(flag_var)
      if(${using_md} MATCHES "True")
        message(WARNING "Petsc was built with /MT, but /MD is currently set.
 See http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F")
      endif(${using_md} MATCHES "True")
    endif (petsc_cc_flags MATCHES "-MT")
  endif (WIN32)

  convert_cygwin_path(petsc_lib_dir)
  message (STATUS "petsc_lib_dir ${petsc_lib_dir}")

  macro (PETSC_FIND_LIBRARY suffix name)
    set (PETSC_LIBRARY_${suffix} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE) # Clear any stale value, if we got here, we need to find it again
    if (WIN32)
      set (libname lib${name}) #windows expects "libfoo", linux expects "foo"
    else (WIN32)
      set (libname ${name})
    endif (WIN32)
    find_library (PETSC_LIBRARY_${suffix} NAMES ${libname} HINTS ${petsc_lib_dir} NO_DEFAULT_PATH)
    set (PETSC_LIBRARIES_${suffix} "${PETSC_LIBRARY_${suffix}}")
    mark_as_advanced (PETSC_LIBRARY_${suffix})
  endmacro (PETSC_FIND_LIBRARY suffix name)

  # Look for petscvec first, if it doesn't exist, we must be using single-library
  petsc_find_library (VEC petscvec)
  if (PETSC_LIBRARY_VEC)
    petsc_find_library (SYS  "petscsys;petsc") # libpetscsys is called libpetsc prior to 3.1 (when single-library was introduced)
    petsc_find_library (MAT  petscmat)
    petsc_find_library (DM   petscdm)
    petsc_find_library (KSP  petscksp)
    petsc_find_library (SNES petscsnes)
    petsc_find_library (TS   petscts)
    macro (PETSC_JOIN libs deps)
      list (APPEND PETSC_LIBRARIES_${libs} ${PETSC_LIBRARIES_${deps}})
    endmacro (PETSC_JOIN libs deps)
    petsc_join (VEC  SYS)
    petsc_join (MAT  VEC)
    petsc_join (DM   MAT)
    petsc_join (KSP  DM)
    petsc_join (SNES KSP)
    petsc_join (TS   SNES)
    petsc_join (ALL  TS)
  else ()
    set (PETSC_LIBRARY_VEC "NOTFOUND" CACHE INTERNAL "Cleared" FORCE) # There is no libpetscvec
    petsc_find_library (SINGLE petsc)
    foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
      set (PETSC_LIBRARIES_${pkg} "${PETSC_LIBRARY_SINGLE}")
    endforeach ()
  endif ()
  if (PETSC_LIBRARY_TS)
    message (STATUS "Recognized Petsc install with separate libraries for each package")
  else ()
    message (STATUS "Recognized Petsc install with single library for all packages")
  endif ()

  include(Check${PETSC_LANGUAGE_BINDINGS}SourceRuns)
  macro (PETSC_TEST_RUNS includes libraries runs)
    if(${PETSC_LANGUAGE_BINDINGS} STREQUAL "C")
      set(_PETSC_ERR_FUNC "CHKERRQ(ierr)")
    elseif(${PETSC_LANGUAGE_BINDINGS} STREQUAL "CXX")
      set(_PETSC_ERR_FUNC "CHKERRXX(ierr)")
    endif()
    if (PETSC_VERSION VERSION_GREATER 3.1)
      set (_PETSC_TSDestroy "TSDestroy(&ts)")
    else ()
      set (_PETSC_TSDestroy "TSDestroy(ts)")
    endif ()

    set(_PETSC_TEST_SOURCE "
static const char help[] = \"Petsc test program.\";
#include <petscts.h>
int main(int argc,char *argv[]) {
  PetscErrorCode ierr;
  TS ts;

  ierr = PetscInitialize(&argc,&argv,0,help);${_PETSC_ERR_FUNC};
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);${_PETSC_ERR_FUNC};
  ierr = TSSetFromOptions(ts);${_PETSC_ERR_FUNC};
  ierr = ${_PETSC_TSDestroy};${_PETSC_ERR_FUNC};
  ierr = PetscFinalize();${_PETSC_ERR_FUNC};
  return 0;
}
")
    multipass_source_runs ("${includes}" "${libraries}" "${_PETSC_TEST_SOURCE}" ${runs} "${PETSC_LANGUAGE_BINDINGS}")
    if (${${runs}})
      set (PETSC_EXECUTABLE_RUNS "YES" CACHE BOOL
        "Can the system successfully run a Petsc executable?  This variable can be manually set to \"YES\" to force CMake to accept a given Petsc configuration, but this will almost always result in a broken build.  If you change PETSC_DIR, PETSC_ARCH, or PETSC_CURRENT you would have to reset this variable." FORCE)
    endif (${${runs}})
  endmacro (PETSC_TEST_RUNS)


  find_path (PETSC_INCLUDE_DIR petscts.h HINTS "${PETSC_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_path (PETSC_INCLUDE_CONF petscconf.h HINTS "${PETSC_DIR}" PATH_SUFFIXES "${PETSC_ARCH}/include" "bmake/${PETSC_ARCH}" NO_DEFAULT_PATH)
  mark_as_advanced (PETSC_INCLUDE_DIR PETSC_INCLUDE_CONF)
  set (petsc_includes_minimal ${PETSC_INCLUDE_CONF} ${PETSC_INCLUDE_DIR})

  petsc_test_runs ("${petsc_includes_minimal}" "${PETSC_LIBRARIES_TS}" petsc_works_minimal)
  if (petsc_works_minimal)
    message (STATUS "Minimal Petsc includes and libraries work.  This probably means we are building with shared libs.")
    set (petsc_includes_needed "${petsc_includes_minimal}")
  else (petsc_works_minimal)     # Minimal includes fail, see if just adding full includes fixes it
    petsc_test_runs ("${petsc_includes_all}" "${PETSC_LIBRARIES_TS}" petsc_works_allincludes)
    if (petsc_works_allincludes) # It does, we just need all the includes (
      message (STATUS "Petsc requires extra include paths, but links correctly with only interface libraries.  This is an unexpected configuration (but it seems to work fine).")
      set (petsc_includes_needed ${petsc_includes_all})
    else (petsc_works_allincludes) # We are going to need to link the external libs explicitly
      resolve_libraries (petsc_libraries_external "${petsc_libs_external}")
      foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
        list (APPEND PETSC_LIBRARIES_${pkg}  ${petsc_libraries_external})
      endforeach (pkg)
      petsc_test_runs ("${petsc_includes_minimal}" "${PETSC_LIBRARIES_TS}" petsc_works_alllibraries)
      if (petsc_works_alllibraries)
         message (STATUS "Petsc only need minimal includes, but requires explicit linking to all dependencies.  This is expected when Petsc is built with static libraries.")
        set (petsc_includes_needed ${petsc_includes_minimal})
      else (petsc_works_alllibraries)
        # It looks like we really need everything, should have listened to Matt
        set (petsc_includes_needed ${petsc_includes_all})
        petsc_test_runs ("${petsc_includes_all}" "${PETSC_LIBRARIES_TS}" petsc_works_all)
        if (petsc_works_all) # We fail anyways
          message (STATUS "Petsc requires extra include paths and explicit linking to all dependencies.  This probably means you have static libraries and something unexpected in Petsc headers.")
        else (petsc_works_all) # We fail anyways
          message (STATUS "Petsc could not be used, maybe the install is broken.")
        endif (petsc_works_all)
      endif (petsc_works_alllibraries)
    endif (petsc_works_allincludes)
  endif (petsc_works_minimal)

  # We do an out-of-source build so __FILE__ will be an absolute path, hence __INSDIR__ is superfluous
  if (${PETSC_VERSION} VERSION_LESS 3.1)
    set (PETSC_DEFINITIONS "-D__SDIR__=\"\"" CACHE STRING "Petsc definitions" FORCE)
  else ()
    set (PETSC_DEFINITIONS "-D__INSDIR__=" CACHE STRING "Petsc definitions" FORCE)
  endif ()
  # Sometimes this can be used to assist FindMPI.cmake
  set (PETSC_MPIEXEC ${petsc_mpiexec} CACHE FILEPATH "Executable for running Petsc MPI programs" FORCE)
  set (PETSC_INCLUDES ${petsc_includes_needed} CACHE STRING "Petsc include path" FORCE)
  set (PETSC_LIBRARIES ${PETSC_LIBRARIES_ALL} CACHE STRING "Petsc libraries" FORCE)
  set (PETSC_COMPILER ${petsc_cc} CACHE FILEPATH "Petsc compiler" FORCE)
  # Note that we have forced values for all these choices.  If you
  # change these, you are telling the system to trust you that they
  # work.  It is likely that you will end up with a broken build.
  mark_as_advanced (PETSC_INCLUDES PETSC_LIBRARIES PETSC_COMPILER PETSC_DEFINITIONS PETSC_MPIEXEC PETSC_EXECUTABLE_RUNS)
endif ()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Petsc
  "Petsc could not be found.  Be sure to set PETSC_DIR and PETSC_ARCH."
  PETSC_INCLUDES PETSC_LIBRARIES PETSC_EXECUTABLE_RUNS)
