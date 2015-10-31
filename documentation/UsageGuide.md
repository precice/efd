# Usage Guide

The project has up to two executable target:

  * `Fluid` --- fluid solver;
  * `RigidBody` --- rigid body solver, `RIGID_BODY` option of CMake must be
    enabled.

The fluid solver must be executed throughout MPI runner (e.g. `mpirun`) with
necessary arguments passed (e.g. `-n 1`).

## Fluid Solver Arguments

  * `--debug,d` --- enable extra debug output from the fluid solver.
  * `--help, -h` --- produce a help message and leave the program.
  * `--no-immersed-boundary, -n` --- disable immersed-boundary computation, i.e. disable PreCICE usage.
  * `--output-directory, -o` --- output directory path, default value: 'output'.
  * `--petsc, -e` --- PETSc configuration path, default value: 'petsc-config.xml'.
  * `--precice, -p` --- PreCICE configuration path, default value: 'precice-config.xml'.
  * `--simulation, -s` --- fluid solver configuration path, default value: 'fluid-config.xml'.

## Fluid Solver XML configuration

    <?xml version="1.0" encoding="utf-8"?>
    <scenario
      <!-- Reynolds number -->
      Re="100"
      <!-- multiplier of the diffusive term in the Navier-Stokes equation
           - could be skiped then 1/Re is used -->
      diffusion-multiplier="0.001"
      <!-- maximum time for the simulation -->
      timeLimit="10.0"
      <!-- maximum iteration number
           - 0 means no limit -->
      iterationLimit="0"
      <!-- time interval for the simulation output
           - 0.0 means every iteration -->
      plotInterval="0.0"
      <!-- safety parameter for the adaptive computation of iteration time step 
           - from 0.0 to 1.0 -->
      tau="0.7"
      <!-- parameter lies between 0 and 1
           - for 0 we recover the central difference discretization,
           and for 1, a pure donor-cell scheme results -->
      gamma="0.5"
      <!-- geometrical width of the domain
           - 2D or 3D vector of floating point values -->
      width="10.0 4.0"
      <!-- cell numbers in the domain
           - 2D or 3D vector of integer values -->
      size="400 160"
      <!-- numbers of processes
           - 2D or 3D vector of integer values -->
      parallelizationSize="1 1"
      <!-- body forces (e.g. gravity force) acting throughout the bulk of the fluid
           - 2D or 3D vector of floating point values -->
      environment="0 0"
      <!-- name prefix for output files -->
      filename="Channel"
      <!-- type of the floating point number
           - float, double, or long double -->
      scalar="double"
      <!-- type of the discretization scheme
           - "Improved Fractional Step Finite Difference" or "Simple Fractional Step Finite Difference" -->
      solver="Improved Fractional Step Finite Difference"
      <!-- type of the output format
           - "Xdmf" or "Vtk" -->
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
      <immersed-boundary
        <!-- iteration from which immersed boundary starts to work
             - This parameter allows the solver to stabalize the value fields
             before applying the immersed-boundary scheme.
             - if omitted it is 0 -->
        start-iteration="10"
        <!-- type of the immersed-boundary algorithm
             - "Precice-Based" or "Rbf-Based" are allowed -->
        type="Precice-Based"
        <!-- enable implicit velicity prediction, otherwise explicit is used
             - if omitted it is 'off' -->
        full-prediction="on"
        <!-- enable moving body,
             the data about its location in the domain is computed begore every iteration
             otherwise body is assumed not moving and
             the data about its location in the domain is computed once
             - if omitted it is 'off' -->
        developing-structure="on"
        <!-- enables computations of the forces acting on the body and send them
             to the structure solver through the PreCICE API
             - This paramter is not necessary if the body is moving in some
             predescribed way, i.e. not related to the fluid flow.
             - if omitted it is 'off' -->
        coupling="on"
        <!-- name of the mesh of the structure body, PreCICE config related -->
        structure-mesh-name="SOLIDZ_Mesh"
        <!-- name of the mesh of the immersed boundary, PreCICE config related -->
        ib-structure-mesh-name="SOLIDZ_Mesh2"
        <!-- name of the data of the acting forces, PreCICE config related -->
        coupling-forces-name="Forces"
        <!-- name of the data of the structure displacements, PreCICE config related -->
        structure-dispacements-name="DisplacementDeltas"
        <!-- name of the data of the immersed-boundary forces, PreCICE config related -->
        immersed-boundary-forces-name="IbForces"
        <!-- file path to the precice config
             - if the path is relative, then it is evaluated from the current XML config path
             - could be omitted then the path from the 'precice' argument of
             the fluid solver is used -->
        precice-configuration-path="../Precice/precice-config.xml"
        <!-- number of cells in outer layers from the body boundary that is used
             to represent the body in the discretized fluid domain -->
        outerLayerSize="0"
        <!-- number of cells in inner layers from the body boundary that is used
             to represent the body in the discretized fluid domain -->
        innerLayerSize="1"/>
    </scenario>
