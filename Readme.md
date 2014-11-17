# Fluid-Structure Interaction (FSI) Simulation

# Functional characteristics

1. Simulate FSI scenarios.
    1. Provide 2D, 3D simulations.
    2. The simulation is fixed-grid.
    3. The fluid solver is finite-difference discretization of the
       incompressible Navier-Stokes equations.
    4. Use different structure solvers:
        1. The rigid-body solver, Structure0815.
        2. *(under consideration)* OpenFOAM solver;
        3. *(under consideration)* COMSOL solver;
        4. *(under consideration)* Alya SOLIDZ solver;
        5. *(under consideration)* CARAT solver.
    5. Continuously resolve boundaries, employ a second-order method.
2. Produce Paraview-readable simulation results.

# Non-functional characteristics

1. Have capability to run in massive parallelization.
2. Have capability to run in sequence.
3. Run on the Linux platform.

# Detailed documentation

- For user documentation, go to the 'documentation/user' directory of the
  project.
- For development documentation, go to the 'documentation/development'
  directory of the project.

