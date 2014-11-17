/*
 * TurbulencePrandtlViscosityStencil.cpp
 *
 *  Created on: Dec 13, 2013
 *      Author: ckow
 */

#include "TurbulencePrandtlViscosityStencil.h"
#include "../stencils/TurbulenceStencilFunctions.h"

TurbulencePrandtlViscosityStencil::TurbulencePrandtlViscosityStencil(
		Parameters & parameters) :
		ComputeLocalViscosityStencil(parameters), _kappa(parameters.blm.kappa), _delta99(
				parameters.blm.delta99) {

	std::cout << "Turbulence model with kappa and delta 99 being\t" << _kappa << '\t' << _delta99 << std::endl;

}

TurbulencePrandtlViscosityStencil::~TurbulencePrandtlViscosityStencil() {}


void TurbulencePrandtlViscosityStencil::apply(TurbulentFlowField & flowField,
int i, int j) {


  loadLocalVelocity2D(flowField, _localVelocity, i, j);
  loadLocalMeshsize2D(_parameters,_localMeshsize,i, j);

  // determine mixing length
  const FLOAT mixingLength = fmin(flowField.getDistanceToWallField().getScalar(i, j) * _kappa,0.09 * _delta99);

  // compute turbulent viscosity
  const FLOAT buffer = mixingLength * mixingLength * sqrt(2.0 * StressTensorSum2D(_localVelocity, _localMeshsize));

  // set turbulent viscosity
  flowField.getViscosityField().getScalar(i, j) = buffer;
}

void TurbulencePrandtlViscosityStencil::apply(TurbulentFlowField & flowField, int i, int j, int k) {
  loadLocalVelocity3D(flowField, _localVelocity, i, j, k);
  loadLocalMeshsize3D(_parameters,_localMeshsize,i, j, k);
  const FLOAT mixingLength = fmin(flowField.getDistanceToWallField().getScalar(i, j, k) * _kappa, 0.09 * _delta99);
  const FLOAT buffer = mixingLength * mixingLength * sqrt(2.0 * StressTensorSum3D(_localVelocity, _localMeshsize));
  flowField.getViscosityField().getScalar(i, j, k) = buffer;
}

