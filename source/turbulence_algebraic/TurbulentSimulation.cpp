#include "TurbulentSimulation.hpp"
//
//
using FsiSimulation::TurbulentSimulation;

// initialise the Simulation-object as usual AND initialise the reference to the
// turbulent flow field data
// as well as the new stencils and iterators
TurbulentSimulation::
TurbulentSimulation(Parameters&         parameters,
                    TurbulentFlowField& flowField)
  : Simulation(parameters, flowField), _stencilFactory(parameters),
    _turbulentFlowField(flowField),
    _computeLocalViscosityStencil(_stencilFactory.getLocalViscosityStencil()),
    _computeLocalViscosityIterator(flowField, parameters,
                                   *_computeLocalViscosityStencil),
    _maxTurbViscStencil(parameters),
    _maxTurbViscIterator(flowField, parameters, _maxTurbViscStencil),
    _turbulentFGHStencil(parameters),
    _turbulentFGHIterator(flowField, parameters, _turbulentFGHStencil),
    _turbulentParallelManager(flowField, parameters)
{}

TurbulentSimulation::
~TurbulentSimulation() {
  // delete local viscosity stencil
  if (_computeLocalViscosityStencil != NULL) {
    delete _computeLocalViscosityStencil;
    _computeLocalViscosityStencil = NULL;
  }
}

void
TurbulentSimulation::
initializeTurbulenceFlowfield() {
  // compute limits of step
  const FLOAT xLimit = _parameters.bfStep.xRatio * _parameters.geometry.lengthX;
  const FLOAT yLimit = _parameters.bfStep.yRatio * _parameters.geometry.lengthY;

  // TODO initialize by iterator and not directly
  if (_parameters.geometry.dim == 2) {
    for (int j = 2; j <= _parameters.parallel.localSize[1] + 1; j++) {
      for (int i = 2; i <= _parameters.parallel.localSize[0] + 1; i++) {
        // if this is a fluid cell...
        if ((_turbulentFlowField.getFlags().getValue(i, j) & OBSTACLE_SELF) ==
            0) {
          // determine offset from step
          FLOAT       stepOffset = 0.0;
          const FLOAT dx         = _parameters.meshsize->getDx(i, j);
          const FLOAT dy         = _parameters.meshsize->getDy(i, j);
          const FLOAT posX       = _parameters.meshsize->getPosX(i, j);
          const FLOAT posY       = _parameters.meshsize->getPosY(i, j);

          if (posX + 0.5 * dx < xLimit) {
            stepOffset = yLimit;
          }

          // determine minimal distance from either bottom or top
          const FLOAT distance = fmin(
            // distance from lower boundary: TODO CHECK
            posY + 0.5 * dy - stepOffset,
            // distance from upper boundary: TODO CHECK
            _parameters.geometry.lengthY - (posY + 0.5 * dy)
            );

          // set distance
          _turbulentFlowField.getDistanceToWallField().getScalar(i, j) =
            distance;
          // for obstacle cells, set distance to zero
        } else {
          _turbulentFlowField.getDistanceToWallField().getScalar(i, j) = 0.0;
        }
      }
    }
  } else {
    for (int k = 2; k <= _parameters.parallel.localSize[2] + 1; k++) {
      for (int j = 2; j <= _parameters.parallel.localSize[1] + 1; j++) {
        for (int i = 2; i <= _parameters.parallel.localSize[0] + 1; i++) {
          // if this is a fluid cell...
          if ((_turbulentFlowField.getFlags().getValue(i, j, k) &
               OBSTACLE_SELF) == 0) {
            // determine offset from step
            FLOAT       stepOffset = 0.0;
            const FLOAT dx         = _parameters.meshsize->getDx(i, j, k);
            const FLOAT dy         = _parameters.meshsize->getDy(i, j, k);
            const FLOAT dz         = _parameters.meshsize->getDz(i, j, k);
            const FLOAT posX       = _parameters.meshsize->getPosX(i, j, k);
            const FLOAT posY       = _parameters.meshsize->getPosY(i, j, k);
            const FLOAT posZ       = _parameters.meshsize->getPosZ(i, j, k);

            if (posX + 0.5 * dx < xLimit) {
              stepOffset = yLimit;
            }

            // determine minimal distance from either bottom or top
            const FLOAT distanceY = fmin(
              // distance from lower boundary: TODO CHECK
              posY + 0.5 * dy - stepOffset,
              // distance from upper boundary: TODO CHECK
              _parameters.geometry.lengthY - (posY + 0.5 * dy)
              );
            const FLOAT distanceZ = fmin(
              // distance from lower boundary: TODO CHECK
              posZ + 0.5 * dz,
              // distance from upper boundary: TODO CHECK
              _parameters.geometry.lengthZ - (posZ + 0.5 * dz)
              );

            // set distance
            _turbulentFlowField.getDistanceToWallField().getScalar(i, j, k) =
              fmin(distanceY, distanceZ);
            // for obstacle cells, set distance to zero
          } else {
            _turbulentFlowField.getDistanceToWallField().getScalar(i, j, k) =
              0.0;
          }
        }
      }
    }
  }
}

void
TurbulentSimulation::
solveTimestep() {
  // compute and set time step
  setTimeStep();

  // compute turbulent viscosity
  _computeLocalViscosityIterator.iterate();
  // exchange viscosity information between neighbouring processes
  _turbulentParallelManager.communicateViscosity();

  // _turbulentFlowField.getViscosityField().show("viscosity");
  // _turbulentFlowField.getDistanceToWallField().show("distance");
  // _turbulentFlowField.getPressure().show("pressure");

  // compute F,G,H based on shear and turbulent viscosity
  _turbulentFGHIterator.iterate();

  // set global boundary values
  _wallFGHIterator.iterate();

  // compute the right hand side
  _rhsIterator.iterate();

  // solve for pressure
  _solver.solve();

  // communicate pressure -> using _parallelManager should also work in this
  // case
  _turbulentParallelManager.communicatePressure();

  // compute velocity
  _velocityIterator.iterate();

  // communicate velocity -> using _parallelManager should also work in this
  // case
  _turbulentParallelManager.communicateVelocity();

  // _turbulentFlowField.getVelocity().show("vel");

  // Iterate for velocities on the boundary
  _wallVelocityIterator.iterate();
}
