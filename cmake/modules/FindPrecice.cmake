
# Precice root directory {{{
# ============================================================================
if (NOT Precice_DIR)
  set(Precice_DIR $ENV{Precice_DIR})
endif()
# ============================================================================
# }}} Precice root directory

# Include {{{
# ============================================================================
if (NOT Precice_INCLUDE_DIRS)
  set(Precice_INCLUDE_DIRS $ENV{Precice_INCLUDE_DIRS})
endif()

if (Precice_INCLUDE_DIRS)
  find_path(PRECICE_INCLUDE_DIR
    NAMES precice/MeshHandle.hpp
    PATHS ${Precice_INCLUDE_DIRS}
    NO_DEFAULT_PATH
    )
endif()

if (Precice_DIR)
  find_path(PRECICE_INCLUDE_DIR
    NAMES src/precice/MeshHandle.hpp
    PATHS ${Precice_DIR}
    NO_DEFAULT_PATH
    )
endif()

if (NOT PRECICE_INCLUDE_DIR)
  find_path(PRECICE_INCLUDE_DIR
    NAMES precice/MeshHandle.hpp
    PATH_SUFFIXES src
    )
endif()

if (PRECICE_INCLUDE_DIR)
  set(PRECICE_INCLUDE_DIRS ${PRECICE_INCLUDE_DIR})
endif()
# ============================================================================
# }}} Include

# Library {{{
# ============================================================================
if (NOT Precice_LIBRARY_DIRS)
  set(Precice_LIBRARY_DIRS $ENV{Precice_LIBRARY_DIRS})
endif()

if (Precice_LIBRARY_DIRS)
  find_library(PRECICE_LIBRARY
    NAMES precice
    PATHS ${Precice_LIBRARY_DIRS}
    PATH_SUFFIXES build/release
    NO_DEFAULT_PATH
    )
endif()

if (Precice_DIR)
  find_library(PRECICE_LIBRARY
    NAMES precice
    PATHS ${Precice_DIR}
    NO_DEFAULT_PATH
    )
endif()

if (NOT PRECICE_LIBRARY)
  find_library(PRECICE_LIBRARY
    NAMES precice
    PATH_SUFFIXES build/release
    )
endif()

if (PRECICE_LIBRARY)
  set(PRECICE_LIBRARIES ${PRECICE_LIBRARY})
endif()
# ============================================================================
# }}} Library

# Finalize {{{
# ============================================================================
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  "PreCICE"
  REQUIRED_VARS
  PRECICE_INCLUDE_DIRS
  PRECICE_LIBRARIES
  )

mark_as_advanced(
  PRECICE_INCLUDE_DIR
  PRECICE_LIBRARY
  )
# ============================================================================
# }}} Finalize
