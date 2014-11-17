/*
 * TurbFGHStencil.cpp
 *
 *  Created on: Dec 13, 2013
 *      Author: ckow
 */

#include "TurbFGHStencil.h"
#include "../stencils/TurbulenceStencilFunctions.h"

TurbFGHStencil::TurbFGHStencil(const Parameters & parameters) : FieldStencil<TurbulentFlowField>(parameters) {}



void TurbFGHStencil::apply(TurbulentFlowField & flowField, int i, int j) {

	// Load local velocities into the center layer of the local array

	loadLocalVelocity2D(flowField, _localVelocity, i, j);
	loadLocalViscosity2D(flowField, _localViscosity, i, j);
        loadLocalMeshsize2D(_parameters,_localMeshsize, i, j);

	FLOAT* const values = flowField.getFGH().getVector(i, j);

	// Now the localVelocity array should contain lexicographically ordered elements around the
	// given index

	// computeF2D like computation
	values[0] = computeTurbulentF2D(_localVelocity, _localMeshsize, _parameters,
			_localViscosity, _parameters.timestep.dt);
	values[1] = computeTurbulentG2D(_localVelocity, _localMeshsize, _parameters,
			_localViscosity, _parameters.timestep.dt);

}

void TurbFGHStencil::apply(TurbulentFlowField & flowField, int i, int j, int k) {
	// The same as in 2D, with slight modifications

	const int obstacle = flowField.getFlags().getValue(i, j, k);

	FLOAT * const values = flowField.getFGH().getVector(i, j, k);

	if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid

		loadLocalVelocity3D(flowField, _localVelocity, i, j, k);
		loadLocalViscosity3D(flowField, _localViscosity, i, j, k);
                loadLocalMeshsize3D(_parameters,_localMeshsize,i,j,k);

		if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
			values[0] = computeTurbulentF3D(_localVelocity, _localMeshsize, _parameters,
					_localViscosity, _parameters.timestep.dt);
		}
		if ((obstacle & OBSTACLE_TOP) == 0) {
			values[1] = computeTurbulentG3D(_localVelocity, _localMeshsize, _parameters,
					_localViscosity, _parameters.timestep.dt);
		}
		if ((obstacle & OBSTACLE_BACK) == 0) {
			values[2] = computeTurbulentH3D(_localVelocity, _localMeshsize, _parameters,
					_localViscosity, _parameters.timestep.dt);
		}
	}
}

