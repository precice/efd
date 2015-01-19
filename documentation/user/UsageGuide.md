# Usage Guide

After installation of the compiled project in the target directory 'bin'
directory appears.
The executable called 'Simulation'.
It is required to be executed throughout MPI runner (e.g. 'mpirun') with
necessary arguments passed (e.g. '-np 1').

For example,

    mpirun -np 1 Simulation

## Arguments

The notation for the arguments is the following: `<full name>, <short name>`.
The full name requires two dashes (`--`), the short --- one dash (`-`).

  * `help, h` --- produce help message.
  * `output-directory, o` --- set output directory path.

    > **NOTE:**
    > It is not required to create directories for the output directory path, it
    > will be created automatically in recursive fashion.

  * `precice, p` --- set PreCiCe configuration path.
  * `simulation, s` --- set fluid simulation configuration path.
  * `petsc, e` --- set PETSc configuration path.


For example,

    mpirun -np 4 ./.install/Release/bin/Simulation -s ./.install/Release/bin/FluidSimulation/Channel.xml -o ~/MyOutputDirector/Channel2DSimpleParabolicInput

