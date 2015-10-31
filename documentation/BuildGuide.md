# Build Guide

# Prerequisites

* For running the build:
    1. GCC 4.9.0 >=; Clang 3.5 >= (recommended);
    2. CMake: version 2.8.8 >=;
* For compiling:
    1. PreCiCe;
    2. PETSc: version 3.3, with enabled MPI support.
    3. Eigen: version 3.2 >=
    4. Boost: version 1.57 >=; with locale, filesystem, program_options, regex, system, thread 
    5. Python libraries: version 2.7
    6. LibXml2
    7. MPI
    8. HDF5

# Quick Instructions

## Building instructions

Here a few Python wrappers are used to simplify CMake usage.

### Steps

1.  Execute CMake configuration with `configure` command.
    - By default generates Unix Makefiles scripts, could be changed with `-G
      <generator>`.
    - By default uses 'build' directory in the `configure` script directory as a
      build directory.
    - By default configures 'install' directory in the `configure` script
      directory as an install directory, could be changed with `-prefix
      <installation path>`
    - By default configures release build, could be changed with `-c
      <configuration type>`, where configuration types are `Release` or `Debug`.
    - By default release mode configuration does not use aggressive compile optimizations
      in favour to short compile time, to enable more hard compile optimizations pass
      `BUILD_WITH_MAX_OPTIMIZATIONS=ON` to the script.
    - By default debug mode configuration does not use advanced debug compile options
      in favour to short compile time, to enable more thorough debuging pass
      `BUILD_WITH_MAX_DEBUGING_ANALYSIS=ON` to the script.
    - To build rigid body solver which code resides in `tools/RigidBody` pass
      `RIGID_BODY=ON` to the script, by default is disabled.

    For example:

        ./configure

    or

        ./configure -G "Ninja"

    or

        ./configure -G "Ninja"
                    -c Debug                            # Release or Debug
                    -prefix ../install                  # Installation path
                    ANY_CMAKE_FLAG=VALUE                # Pass flags to CMake,
                    ...                                 # will be prefixed by '-D'
                    ...                                 # automatically
                    BUILD_WITH_MAX_DEBUGING_ANALYSIS=ON # Enables additional debuging options
                    RIGID_BODY=ON                       # Enables rigid body solver

2. Execute building with `compile` command.
    - By default, use Unix make tool.

    For example:

        ./compile

    or to use `ninja` as a build tool

        ./compile ninja

     or set the `V` variable to `1`, what, in case of the Unix Makefiles,
     makes output verbose:

        ./compile V=1 # Make Unix Makefiles verbose

     or pass any flags to the `compile` command that must be passed to the low-level
     build system command:

        ./compile -MyFlagToBuildSystem # Any flags will be passed
                  ...                  # to build system command

# Main Instructions

## Building instructions

The common procedure for CMake is used, therefore for more details go to the
CMake user guide.

### Steps

1.  Execute configuration CMake command.
    - By default release mode configuration does not use aggressive compile optimizations
      in favour to short compile time, to enable more hard compile optimizations pass
      `BUILD_WITH_MAX_OPTIMIZATIONS=ON` to the script.
    - By default debug mode configuration does not use advanced debug compile options
      in favour to short compile time, to enable more thorough debuging pass
      `BUILD_WITH_MAX_DEBUGING_ANALYSIS=ON` to the script.
    - To build rigid body solver which code resides in `tools/RigidBody` pass
      `RIGID_BODY=ON` to the script, by default is disabled.

    For example:

        cmake -G "Unix Makefiles"
              -DCMAKE_BUILD_TYPE=Release            # Release or Debug
              -DCMAKE_INSTALL_PREFIX=../install     # Installation Path
              -DBUILD_WITH_MAX_OPTIMIZATIONS=ON     # Enables additional optimizations
              -DRIGID_BODY=ON                       # Enables rigid body solver
              ../

2.  Execute building command.

    For example:

        make

### Dependency lookup

The build system will try to find dependencies by itself, for more information
go to the CMake manual (`find_path` for includes; `find_library` for libraries).

> **NOTE:**
> MPI implementation must be the same that was used to compile PETSc.

> **NOTE:**
> Python libraries must be the same that was used to compile PreCICE.

#### Configuration of the dependency lookup

If the dependency are not visible throughout environment, or properly install in
the system, then use the instructions below before CMake configuration for a
corresponding dependency.

All the dependencies are configured in the same fashion.

For each dependency there are variable names.
Use this variable names in one or both ways to control the search of a
dependency on the system:

+ as environment variables;
+ as CMake configuration variables.

The lookup of the following dependencies could be configured:

*   **PETSc**

    + `PETSC_DIR` --- directory in which PETSc resides.
    + *(optional)* `PETSC_ARCH` --- build architecture.

*   **MPI**

    + `MPI_HOME` --- directory, where MPI resides.

*   **PreCICE**

    + *(optional)* `PRECICE_INCLUDE_DIRS` --- directories, where PreCICE's includes
      resides.
    + *(optional)* `PRECICE_LIBRARY_DIRS` --- directories, where PreCICE's libraries
      resides.
    + *(optional)* `PRECICE_DIR` --- directory, where PreCICE resides.

    Using variables above the order of search is described below, it moves to
    the next step until the check of some conditions succeed.

    Lookup for includes in values of the variables in the order:

    1. `PRECICE_INCLUDE_DIRS`
    2. `PRECICE_DIR`

    Lookup for libraries in values of the variables in the order:

    1. `PRECICE_LIBRARY_DIRS`
    2. `PRECICE_DIR`

*   **Eigen**

    + `EIGEN_INCLUDE_DIRS` --- directories, where Eigen's includes resides.

*   **Boost**

    + `BOOST_ROOT` (or `BOOSTROOT`) --- preferred installation prefix.
    + `BOOST_INCLUDEDIR` --- preferred include directory e.g. '<prefix>/include'.
    + `BOOST_LIBRARYDIR` --- preferred library directory e.g. '<prefix>/lib'.
    + `Boost_NO_SYSTEM_PATHS` --- set to 'ON' to disable searching in locations
      not specified by these hint variables. Default is 'OFF'.
    + `Boost_ADDITIONAL_VERSIONS` --- list of Boost versions not known to this
      module (Boost install locations may contain the version).

*   **Python libraries**

    + `PYTHON_LIBRARY` --- path to the python library.
    + `PYTHON_INCLUDE_DIR` --- path to where Python.h is found.
