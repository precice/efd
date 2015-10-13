# General {{{
# ==============================================================================
cmake_minimum_required(VERSION 2.8.8)
# ------------------------------------------------------------------------------
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif ()
# ------------------------------------------------------------------------------
project("Fluid-Structure Interations Simulation" CXX C)
# ------------------------------------------------------------------------------
set(CMAKE_ERROR_DEPRECATED ON)
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
# General }}}

# Options {{{
# ==============================================================================
option(RIGID_BODY
  "Compiles rigid body solver from 'tools/RigidBody'" OFF)
option(BUILD_WITH_MAX_OPTIMIZATIONS
  "Add additional flags to the compiler to provide maximum optimizations" OFF)
option(BUILD_WITH_MAX_DEBUGING_ANALYSIS
  "Add additional flags to the compiler to provide maximum debuging analysis" OFF)
# ==============================================================================
# Options }}}

# Modules {{{
# ==============================================================================
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)
# ==============================================================================
# Modules }}}

# Packages {{{
# ==============================================================================
find_package(Petsc REQUIRED)
if ("${MPI_CXX_COMPILER}" STREQUAL "")
  set(MPI_CXX_COMPILER ${PETSC_COMPILER})
endif()
find_package(MPI REQUIRED)
find_package(Precice REQUIRED)
find_package(PythonLibs 2.7 REQUIRED)
set(Boost_USE_STATIC_LIBS    OFF)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.57
  COMPONENTS
  locale
  filesystem
  program_options
  regex
  system
  thread # Required for precice
  REQUIRED)
find_package(Eigen REQUIRED)
find_package(LibXml2 REQUIRED)
find_path(Hdf5_INCLUDE_DIR
  hdf5.h
  HINTS $ENV{HDF5_ROOT}
  PATH_SUFFIXES include
)
find_library(Hdf5_LIBRARY
  hdf5
  HINTS $ENV{HDF5_ROOT}
  PATH_SUFFIXES lib
)
# ==============================================================================
# Packages }}}

# Compilers {{{
# ==============================================================================
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(COMPILER_IS_CLANG TRUE)
  set(LINKER_IS_LD      TRUE)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(COMPILER_IS_GCC TRUE)
  set(LINKER_IS_LD    TRUE)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(COMPILER_IS_ICC TRUE)
  set(LINKER_IS_ICC    TRUE)
else()
  message(FATAL "Unsopported toolchain")
endif ()
# ==============================================================================
# Compilers }}}

# {{{ General-purpose definitions
# ==============================================================================
# ==============================================================================
# General-purpose definitions }}}

# General-purpose flags {{{
# ==============================================================================
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if (COMPILER_IS_CLANG)
  # Hide unused compiler flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
elseif (COMPILER_IS_GCC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
elseif (COMPILER_IS_ICC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ip -ipo2 -ipo-jobs4")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -opt-prefetch=2")
endif ()
if (COMPILER_IS_CLANG OR COMPILER_IS_GCC)
  if (BUILD_WITH_MAX_OPTIMIZATIONS)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
  endif()
endif()
if (COMPILER_IS_CLANG OR COMPILER_IS_GCC OR COMPILER_IS_ICC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

  if (BUILD_WITH_MAX_OPTIMIZATIONS)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops")
  else()
    string(REPLACE "-O3" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O0")
  endif()

  if (BUILD_WITH_MAX_DEBUGING_ANALYSIS)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-inline")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=memory")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
endif()
endif ()
# ------------------------------------------------------------------------------
if (LINKER_IS_LD)
  if (BUILD_WITH_MAX_OPTIMIZATIONS)
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,-O5")
  endif()
endif ()
# ==============================================================================
# General-purpose flags }}}

# Targets {{{
# ==============================================================================
file(GLOB_RECURSE UNI_CPP_FILES "${PROJECT_DIR}/3rdParty/Uni/source/*.cpp")
add_library(Uni STATIC ${UNI_CPP_FILES})
target_include_directories(Uni PUBLIC "${PROJECT_DIR}/3rdParty/Uni/include")
target_include_directories(Uni SYSTEM PRIVATE ${Boost_INCLUDE_DIRS} ${Eigen_INCLUDE_DIRS})
target_link_libraries(Uni PRIVATE ${Boost_LIBRARIES})

file(GLOB_RECURSE FLUID_CPP_FILES
  "${PROJECT_DIR}/source/*.cpp"
  )

set(Fluid_Private_INCLUDE_DIRS
  "${PROJECT_DIR}/source"
  )

set(Fluid_INCLUDES
  ${Precice_INCLUDE_DIRS}
  ${PETSC_INCLUDES}
  ${Eigen_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${LIBXML2_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${Hdf5_INCLUDE_DIR}
  )
set(Fluid_DEFINITIONS
  -DMPICH_IGNORE_CXX_SEEK
  -DPETSC_USE_REAL_DOUBLE
  ${PETSC_DEFINITIONS}
  ${LIBXML2_DEFINITIONS}
  ${HDF5_DEFINITIONS}
  )
set(Fluid_COMPILE_FLAGS
  ${MPI_CXX_COMPILE_FLAGS}
  )
add_executable(Fluid ${FLUID_CPP_FILES})
target_include_directories(Fluid PRIVATE ${Fluid_Private_INCLUDE_DIRS})
target_include_directories(Fluid SYSTEM PRIVATE ${Fluid_INCLUDES})
target_link_libraries(Fluid PRIVATE
  Uni
  ${Precice_LIBRARIES}
  ${PETSC_LIBRARIES}
  ${Boost_LIBRARIES}
  ${LIBXML2_LIBRARIES}
  ${PYTHON_LIBRARIES}
  ${Hdf5_LIBRARY}
  )

if (RIGID_BODY)
  file(GLOB_RECURSE STRUCTURE_CPP_FILES "${PROJECT_DIR}/tools/RigidBody/source/*.cpp")
  add_executable(Structure ${STRUCTURE_CPP_FILES})
  target_include_directories(Structure PRIVATE "${PROJECT_DIR}/tools/RigidBody/source")
  target_include_directories(Structure SYSTEM PRIVATE
    ${Precice_INCLUDE_DIRS}
    ${Eigen_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${LIBXML2_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_PATH}
    )
  target_link_libraries(Structure PRIVATE
  Uni
    ${Precice_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${Boost_LIBRARIES}
    ${LIBXML2_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    )
endif()
# ==============================================================================
# Targets }}}

# Targets' definitions {{{
# ==============================================================================
target_compile_definitions(Fluid PRIVATE ${Fluid_DEFINITIONS})
if (RIGID_BODY)
  target_compile_definitions(Structure PRIVATE ${LIBXML2_DEFINITIONS})
endif()
# ==============================================================================
# Targets' definitions }}}

# Targets' flags {{{
# ==============================================================================
if(MPI_CXX_COMPILE_FLAGS)
  set_target_properties(Fluid
    PROPERTIES COMPILE_FLAGS ${Fluid_COMPILE_FLAGS})
  if (RIGID_BODY)
    set_target_properties(Structure
      PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif()
endif()

if(MPI_CXX_LINK_FLAGS)
  set_target_properties(Fluid
    PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  if (RIGID_BODY)
    set_target_properties(Structure
      PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif()
endif()
# ==============================================================================
# Targets' flags }}}

# Installation {{{
# ==============================================================================
install(
  TARGETS Fluid
  ARCHIVE DESTINATION ${INSTALL_LIBRARY_DIR}
  LIBRARY DESTINATION ${INSTALL_LIBRARY_DIR}
  RUNTIME DESTINATION ${INSTALL_BINARY_DIR}
  )

if (RIGID_BODY)
  install(
    TARGETS Structure
    ARCHIVE DESTINATION ${INSTALL_LIBRARY_DIR}
    LIBRARY DESTINATION ${INSTALL_LIBRARY_DIR}
    RUNTIME DESTINATION ${INSTALL_BINARY_DIR}
    )
endif()
# ==============================================================================
# Installation }}}
