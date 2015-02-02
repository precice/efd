# Precice root directory {{{
# ============================================================================
if (NOT PRECICE_DIR)
  set(PRECICE_DIR $ENV{PRECICE_DIR})
endif()
# ============================================================================
# }}} Precice root directory

# Include {{{
# ============================================================================
if (NOT PRECICE_INCLUDE_DIRS)
  set(PRECICE_INCLUDE_DIRS $ENV{PRECICE_INCLUDE_DIRS})
endif()

if (PRECICE_INCLUDE_DIRS)
  find_path(Precice_INCLUDE_DIR
    NAMES precice/MeshHandle.hpp
    PATHS ${PRECICE_INCLUDE_DIRS}
    NO_DEFAULT_PATH
    )
endif()

if (PRECICE_DIR)
  find_path(Precice_INCLUDE_DIR
    NAMES precice/MeshHandle.hpp
    PATHS ${PRECICE_DIR}
    PATH_SUFFIXES src
    NO_DEFAULT_PATH
    )
endif()

if (NOT Precice_INCLUDE_DIR)
  find_path(Precice_INCLUDE_DIR
    NAMES precice/MeshHandle.hpp
    PATH_SUFFIXES src
    )
endif()

if (Precice_INCLUDE_DIR)
  set(Precice_INCLUDE_DIRS ${Precice_INCLUDE_DIR})
endif()
# ============================================================================
# }}} Include

# Library {{{
# ============================================================================
if (NOT PRECICE_LIBRARY_DIRS)
  set(PRECICE_LIBRARY_DIRS $ENV{PRECICE_LIBRARY_DIRS})
endif()

if (PRECICE_LIBRARY_DIRS)
  find_library(Precice_LIBRARY
    NAMES precice
    PATHS ${PRECICE_LIBRARY_DIRS}
    NO_DEFAULT_PATH
    )
endif()

set(_Precice_build_type "release")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(_Precice_build_type "debug")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

if (PRECICE_DIR)
  find_library(Precice_LIBRARY
    NAMES precice
    PATHS ${PRECICE_DIR}
    PATH_SUFFIXES "build/${_Precice_build_type}"
    NO_DEFAULT_PATH
    )
endif()

if (NOT Precice_LIBRARY)
  find_library(Precice_LIBRARY
    NAMES precice
    PATH_SUFFIXES "build/${_Precice_build_type}"
    )
endif()

if (Precice_LIBRARY)
  set(Precice_LIBRARIES ${Precice_LIBRARY})
endif()
# ============================================================================
# }}} Library

# Finalize {{{
# ============================================================================
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  "PreCICE"
  REQUIRED_VARS
  Precice_INCLUDE_DIRS
  Precice_LIBRARIES
  )

mark_as_advanced(
  Precice_INCLUDE_DIR
  Precice_LIBRARY
  )
# ============================================================================
# }}} Finalize
