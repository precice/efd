# General {{{
# ==============================================================================
cmake_minimum_required(VERSION 2.8.8)
# ------------------------------------------------------------------------------
#set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)

#if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
#  message(SEND_ERROR "In-source builds are not allowed.")
#endif ()
# ------------------------------------------------------------------------------
project("Fsi Simulation" C CXX)
# ------------------------------------------------------------------------------
set(CMAKE_ERROR_DEPRECATED ON)
# ------------------------------------------------------------------------------
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE   ON)
# ------------------------------------------------------------------------------
if (WIN32)
  set(CMAKE_SHARED_LIBRARY_PREFIX "")
endif ()
# ------------------------------------------------------------------------------
get_filename_component(PROJECT_DIR "." ABSOLUTE)
set(INSTALL_BINARY_DIR  bin)
set(INSTALL_INCLUDE_DIR include)
set(INSTALL_LIBRARY_DIR lib)
if (WIN32)
  set(INSTALL_CMAKE_DIR cmake)
else()
  set(INSTALL_CMAKE_DIR lib/cmake)
endif ()
# ==============================================================================
# }}} General

# Modules {{{
# ==============================================================================
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)
# ==============================================================================
# }}} Modules

# Packages {{{
# ==============================================================================
#find_package(Petsc COMPONENTS CXX REQUIRED)
find_package(Petsc REQUIRED)
set(MPI_CXX_COMPILER ${PETSC_COMPILER})
find_package(MPI REQUIRED)
find_package(Precice REQUIRED)
find_package(PythonLibs 2.7 REQUIRED)
set(Boost_USE_STATIC_LIBS    OFF)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.57 # Versions 1.54, 1.55 have failed
  COMPONENTS
  locale
  filesystem
  program_options
  regex
  system
  thread # Required for precice
  REQUIRED)
find_package(Precice REQUIRED)
find_package(Eigen REQUIRED)
find_package(Uni REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(XDMF REQUIRED)
find_library(XDMF_LIBRARIES
  NAMES Xdmf
  PATHS ${XDMF_LIBRARY_DIRS}
  NO_DEFAULT_PATH
  )
# ==============================================================================
# }}} Packages

# Compilers {{{
# ==============================================================================
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(COMPILER_IS_G++ TRUE)
  set(LINKER_IS_LD    TRUE)
else()
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(COMPILER_IS_G++ TRUE)
  set(LINKER_IS_LD    TRUE)
else()
  set(COMPILER_IS_G++ TRUE)
  set(LINKER_IS_LD    TRUE)
endif ()
endif ()
# ==============================================================================
# }}} Compilers

# {{{ General-purpose definitions
# ==============================================================================
# ==============================================================================
# }}} General-purpose definitions

# General-purpose flags {{{
# ==============================================================================
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if (COMPILER_IS_G++)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")
 # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fomit-frame-pointer")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops")

  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -fno-inline")
  #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
endif ()
# ------------------------------------------------------------------------------
if (LINKER_IS_LD)
  if (WIN32)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  endif ()

  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,-O5")
  # ----------------------------------------------------------------------------
  if (WIN32)
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -static-libgcc")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -static-libstdc++")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -static")
  endif ()

  set(CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS_RELEASE} -Wl,-O5")
  # ----------------------------------------------------------------------------
  if (WIN32)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static")
  endif ()

  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -Wl,-O5")
  # ----------------------------------------------------------------------------
  set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS_RELEASE} -Wl,-O5")
endif ()
# ==============================================================================
# }}} General-purpose flags

# Targets {{{
# ==============================================================================
file(GLOB_RECURSE CPP_FILES
  "${PROJECT_DIR}/source/**/*.cpp"
  "${PROJECT_DIR}/source/**/*.hpp"
  )

add_executable(Simulation ${CPP_FILES})

target_include_directories(Simulation PRIVATE
  "${PROJECT_DIR}/source"
  )

target_include_directories(Simulation SYSTEM PRIVATE
  ${Precice_INCLUDE_DIRS}
  ${PETSC_INCLUDES}
  ${Eigen_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${LIBXML2_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${XDMF_HDF5_INCLUDE}
  ${XDMF_LIBXML2_INCLUDE}
  ${XDMF_ZLIB_INCLUDE}
  )
target_link_libraries(Simulation PRIVATE
  Uni
  ${Precice_LIBRARIES}
  ${PETSC_LIBRARIES}
  ${Boost_LIBRARIES}
  ${LIBXML2_LIBRARIES}
  ${PYTHON_LIBRARIES}
  ${XDMF_LIBRARIES}
  ${XDMF_LIBXML2_LIBRARIES}
  ${XDMF_HDF5_LIBRARIES}
  ${XDMF_ZLIB_LIBRARIES}
  )
# ==============================================================================
# }}} Targets

# Targets' definitions {{{
# ==============================================================================
target_compile_definitions(Simulation PRIVATE
  ${PETSC_DEFINITIONS}
  ${LIBXML2_DEFINITIONS})
# ==============================================================================
# }}} Targets' definitions

# Targets' flags {{{
# ==============================================================================
if(MPI_CXX_COMPILE_FLAGS)
  set_target_properties(Simulation PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_CXX_LINK_FLAGS)
  set_target_properties(Simulation PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()
# ==============================================================================
# }}} Targets' flags

# {{{ Installation
# ==============================================================================
install(
  TARGETS Simulation
  ARCHIVE DESTINATION ${INSTALL_LIBRARY_DIR}
  LIBRARY DESTINATION ${INSTALL_LIBRARY_DIR}
  RUNTIME DESTINATION ${INSTALL_BINARY_DIR}
  )

install(DIRECTORY "${PROJECT_DIR}/configuration/FluidPetsc"
  DESTINATION "${INSTALL_BINARY_DIR}"
  FILES_MATCHING REGEX "(.*\\.conf)$"
  )

install(DIRECTORY "${PROJECT_DIR}/configuration/Precice"
  DESTINATION "${INSTALL_BINARY_DIR}"
  FILES_MATCHING REGEX "(.*\\.(xml|py))$"
  )


install(DIRECTORY "${PROJECT_DIR}/configuration/FluidSimulation"
  DESTINATION "${INSTALL_BINARY_DIR}"
  FILES_MATCHING REGEX "(.*\\.xml)$"
  )
# ==============================================================================
# }}} Installation
