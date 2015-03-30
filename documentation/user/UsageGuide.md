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


Known issues:

*   Parallel and sequential versions give different results and different
    number of iterations to converge. From [PETSc FAQ][1]

    > **How come when I run the same program on the same number of processes I get
    a "different" answer?**
    >
    > Inner products and norms in PETSc are computed using the MPI_Allreduce() command.In different runs the order at which values arrive at a given process (via MPI) can be in a different order, thus the order in which some floating point arithmetic operations are performed will be different. Since floating point arithmetic arithmetic is not associative, the computed quantity may be (slightly) different. Over a run the many slight differences in the inner products and norms will effect all the computed results. It is important to realize that none of the computed answers are any less right or wrong (in fact the sequential computation is no more right then the parallel ones), they are all equally valid. 
    >
    > The discussion above assumes that the exact same algorithm is being used on the different number of processes. When the algorithm is different for the different number of processes (almost all preconditioner algorithms except Jacobi are different for different number of processes) then one expects to see (and does) a greater difference in results for different numbers of processes. In some cases (for example block Jacobi preconditioner) it may be that the algorithm works for some number of processes and does not work for others.

    > **How come when I run the same linear solver on a different number of
    processes it takes a different number of iterations?**
    >
    > The convergence of many of the preconditioners in PETSc including the the default parallel preconditioner block Jacobi depends on the number of processes. The more processes the (slightly) slower convergence it has. This is the nature of iterative solvers, the more parallelism means the more "older" information is used in the solution process hence slower convergence

[1]: http://www.mcs.anl.gov/petsc/documentation/faq.html
