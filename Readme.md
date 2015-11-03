A Fixed-Grid Flow Solver for Fluid-Structure Interaction with The Coupling Library PreCICE
==========================================================================================

Description
-----------

--------------------------------------------------------------------------------

This project is a C++ implementation of the fluid motion numerical simulation
of the incompressible [Navier-Stokes equations][NS] on a [Eulerian
grid][Eulerian] for modeling of [fluid-structure interaction][FSI] scientific
phenomena with [the partitioned approach][Partition] using the coupling library
[PreCICE][].

Features
--------

--------------------------------------------------------------------------------

* The fluid motion is modeled by the incompressible Navier-Stokes equations.
* For the temporal descritization the projection schemes are used.
* A finite-differnce on a Cartesian grid is used for the spatial descritization.
* Two- and three-dimensional scenarios are available
* The application produces Paraview-readable simulation results in [HDF5][] and [VTK][]
  formats
* The solver can run in series and in parallel using a domain decomposition approach.

--------------------------------------------------------------------------------

### Flow Around a Disk

![Flow Around a Disk #1](https://dl.dropbox.com/s/6zjqwbyvjcv2qe0/u.png)

![Flow Around a Disk #2](https://dl.dropbox.com/s/f4r6nqtf5kv7gxf/v.png)

![Flow Around a Disk #3](https://dl.dropbox.com/s/s8rs3ifh61o3q86/pressure.png)

![Flow Around a Disk #4](https://dl.dropbox.com/s/rfw7fw0v6x1axiq/vorticity.png)

--------------------------------------------------------------------------------

### Predescribed Motion of a Disk

![Predescribed Motion of a Disk #1](https://dl.dropbox.com/s/szmc3udx276zdrx/u1.png)

![Predescribed Motion of a Disk #2](https://dl.dropbox.com/s/rb2w0jfp76pbgtr/v1.png)

![Predescribed Motion of a Disk #3](https://dl.dropbox.com/s/279fsz39c76piuq/pressure1.png)

![Predescribed Motion of a Disk #4](https://dl.dropbox.com/s/7c287hyslez3mwt/vorticity1.png)

![Predescribed Motion of a Disk #5](https://dl.dropbox.com/s/o6jh3nept0upl6z/u2.png)

![Predescribed Motion of a Disk #6](https://dl.dropbox.com/s/vihjac3xtzbvhuh/v2.png)

![Predescribed Motion of a Disk #7](https://dl.dropbox.com/s/8lajz7ybjqjvwjz/pressure2.png)

![Predescribed Motion of a Disk #8](https://dl.dropbox.com/s/0zlxwokuj7x5qnm/vorticity2.png)

--------------------------------------------------------------------------------

### Free Motion of A Disk

![Free Motion of a Disk #1](https://dl.dropbox.com/s/ia3lsv9v64to03w/u.png)

![Free Motion of a Disk #2](https://dl.dropbox.com/s/bxvkacvbb8d4ns6/v.png)

![Free Motion of a Disk #3](https://dl.dropbox.com/s/a8hvtesu9btc6jg/pressure.png)

![Free Motion of a Disk #4](https://dl.dropbox.com/s/1omieo3isikrv6g/vorticity.png)

--------------------------------------------------------------------------------

### Free Fall

![Bubble Scenario #1](https://dl.dropbox.com/s/zqt8b99hrzd7ggf/Free-Fall.gif)

Build Guide
-----------

--------------------------------------------------------------------------------

### Prerequisites

--------------------------------------------------------------------------------

* For running the build:
    1. GCC 4.9.0 >=; Clang 3.5 >= (recommended);
    2. [CMake][]: version 2.8.8 >=;
* For compiling:
    1. [PreCICE][];
    2. [PETSc][]: version 3.3, with enabled MPI support.
    3. [Eigen][]: version 3.2 >=
    4. [Boost][]: version 1.57 >=; with locale, filesystem, program_options, regex, system, thread 
    5. Python libraries: version 2.7
    6. [LibXml2][]
    7. [MPI][]
    8. [HDF5][]

### Quick Building Instructions

--------------------------------------------------------------------------------

Here a few Python wrappers are used to simplify CMake usage in the following steps.

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

### CMake Building Instructions

--------------------------------------------------------------------------------

The common procedure for CMake is used, therefore for more details go to the
CMake user guide.

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

#### Dependency lookup

--------------------------------------------------------------------------------

The build system will try to find dependencies by itself, for more information
go to the CMake manual (`find_path` for includes; `find_library` for libraries).

> **NOTE:**
> MPI implementation must be the same that was used to compile PETSc.

> **NOTE:**
> Python libraries must be the same that was used to compile PreCICE.

#### Configuration of the dependency lookup

--------------------------------------------------------------------------------

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


Usage Guide
-----------

--------------------------------------------------------------------------------

The project has up to two executable target:

  * `Fluid` --- fluid solver;
  * `RigidBody` --- rigid body solver, `RIGID_BODY` option of CMake must be
    enabled.

The fluid solver must be executed throughout MPI runner (e.g. `mpirun`) with
necessary arguments passed (e.g. `-n 1`).

### Fluid Solver Arguments

--------------------------------------------------------------------------------

  * `--debug,d` --- enable extra debug output from the fluid solver.
  * `--help, -h` --- produce a help message and leave the program.
  * `--no-immersed-boundary, -n` --- disable immersed-boundary computation, i.e. disable PreCICE usage.
  * `--output-directory, -o` --- output directory path, default value: 'output'.
  * `--petsc, -e` --- PETSc configuration path, default value: 'petsc-config.xml'.
  * `--precice, -p` --- PreCICE configuration path, default value: 'precice-config.xml'.
  * `--simulation, -s` --- fluid solver configuration path, default value: 'fluid-config.xml'.

### Fluid Solver XML configuration

--------------------------------------------------------------------------------

    <?xml version="1.0" encoding="utf-8"?>
      <!-- `Re`
           Reynolds number -->
      <!-- `diffusion-multiplier`
           multiplier of the diffusive term in the Navier-Stokes equation
           - could be skiped then 1/Re is used -->
      <!-- `timeLimit`
           maximum time for the simulation -->
      <!-- iterationLimit
           maximum iteration number
           - 0 means no limit -->
      <!-- `plotInterval`
           time interval for the simulation output
           - 0.0 means every iteration -->
      <!-- `tau`
           safety parameter for the adaptive computation of iteration time step 
           - from 0.0 to 1.0 -->
      <!-- `gamma`
           parameter lies between 0 and 1
           - for 0 we recover the central difference discretization,
           and for 1, a pure donor-cell scheme results -->
      <!-- `width`
           geometrical width of the domain
           - 2D or 3D vector of floating point values -->
      <!-- `size`
           cell numbers in the domain
           - 2D or 3D vector of integer values -->
      <!-- `parallelizationSize`
           numbers of processes
           - 2D or 3D vector of integer values -->
      <!-- `environment`
           body forces (e.g. gravity force) acting throughout the bulk of the fluid
           - 2D or 3D vector of floating point values -->
      <!-- `filename`
           name prefix for output files -->
      <!-- `scalar`
           type of the floating point number
           - float, double, or long double -->
      <!-- `solver`
           type of the discretization scheme
           - "Improved Fractional Step Finite Difference" or "Simple Fractional Step Finite Difference" -->
      <!-- `output`
           type of the output format
           - "Xdmf" or "Vtk" -->
    <scenario
      Re="100"
      diffusion-multiplier="0.001"
      timeLimit="10.0"
      iterationLimit="0"
      plotInterval="0.0"
      tau="0.7"
      gamma="0.5"
      width="10.0 4.0"
      size="400 160"
      parallelizationSize="1 1"
      environment="0 0"
      filename="Channel"
      scalar="double"
      solver="Improved Fractional Step Finite Difference"
      output="Xdmf"
      >
      <!-- boundary condition specification -->
      <walls>
        <!-- boundary condition for left, right, bottom, top, back, or front
             wall
             - `type` --- attribute that specifies the condition type
                        "ParabolicInput", "Input", or Output are allowed"
             - `velocity` -- atribute that specifies velocity for "ParabolicInput"
                           or "Input" types
                           2D or 3D vector of floating point values
             -->
        <left type="ParabolicInput" velocity="1.5 0.0" />
        <right type="Output"  velocity="0.0 0.0"/>
        <bottom type="Input" velocity="0.0 0.0" />
        <top type="Input" velocity="0.0 0.0" />
      </walls>
      <!-- settings for immersed boundary problem -->
        <!-- `start-iteration`
             iteration from which immersed boundary starts to work
             - This parameter allows the solver to stabalize the value fields
             before applying the immersed-boundary scheme.
             - if omitted it is 0 -->
        <!-- `type`
             type of the immersed-boundary algorithm
             - "Precice-Based" or "Rbf-Based" are allowed -->
        <!-- `full-prediction`
             enable implicit velicity prediction, otherwise explicit is used
             - if omitted it is 'off' -->
        <!-- `developing-structure`
             enable moving body,
             the data about its location in the domain is computed begore every iteration
             otherwise body is assumed not moving and
             the data about its location in the domain is computed once
             - if omitted it is 'off' -->
        <!-- `coupling`
             enables computations of the forces acting on the body and send them
             to the structure solver through the PreCICE API
             - This paramter is not necessary if the body is moving in some
             predescribed way, i.e. not related to the fluid flow.
             - if omitted it is 'off' -->
        <!-- `structure-mesh-name`
             name of the mesh of the structure body, PreCICE config related -->
        <!-- `ib-structure-mesh-name`
             name of the mesh of the immersed boundary, PreCICE config related -->
        <!-- `coupling-forces-name`
             name of the data of the acting forces, PreCICE config related -->
        <!-- `structure-displacements-name`
             name of the data of the structure displacements, PreCICE config related -->
        <!-- `immersed-boundary-forces-name`
             name of the data of the immersed-boundary forces, PreCICE config related -->
        <!-- `precice-configuration-path`
             file path to the precice config
             - if the path is relative, then it is evaluated from the current XML config path
             - could be omitted then the path from the 'precice' argument of
             the fluid solver is used -->
        <!-- `outerLayerSize`
             number of cells in outer layers from the body boundary that is used
             to represent the body in the discretized fluid domain -->
        <!-- `innerLayerSize`
             number of cells in inner layers from the body boundary that is used
             to represent the body in the discretized fluid domain -->
      <immersed-boundary
        start-iteration="10"
        type="Precice-Based"
        full-prediction="on"
        developing-structure="on"
        coupling="on"
        structure-mesh-name="SOLIDZ_Mesh"
        ib-structure-mesh-name="SOLIDZ_Mesh2"
        coupling-forces-name="Forces"
        structure-dispacements-name="DisplacementDeltas"
        immersed-boundary-forces-name="IbForces"
        precice-configuration-path="../Precice/precice-config.xml"
        outerLayerSize="0"
        innerLayerSize="1"/>
    </scenario>

Copyright
---------

--------------------------------------------------------------------------------

Copyright (C) 2015, Viacheslav Mikerov

License
-------

--------------------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.


[NS]:        https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations
[Eulerian]:  https://en.wikipedia.org/wiki/Lagrangian_and_Eulerian_specification_of_the_flow_field
[FSI]:       https://en.wikipedia.org/wiki/Fluid%E2%80%93structure_interaction
[Partition]: https://en.wikipedia.org/wiki/Fluid%E2%80%93structure_interaction#Analysis

[PreCICE]:   https://github.com/precice/precice
[HDF5]:      https://www.hdfgroup.org/HDF5/
[VTK]:       http://www.paraview.org/Wiki/ParaView/Data_formats#VTK.28Visualization_ToolKit.29_files
[CMake]:     https://cmake.org/
[Eigen]:     http://eigen.tuxfamily.org/index.php?title=Main_Page
[Boost]:     http://www.boost.org/
[LibXML2]:   http://www.xmlsoft.org/
[MPI]:       https://en.wikipedia.org/wiki/Message_Passing_Interface 
[PETSc]:     http://www.mcs.anl.gov/petsc/
