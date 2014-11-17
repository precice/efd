/*
 * BaldwinLomaxTurbulentViscosityStencil.cpp
 *
 *  Created on: Dec 14, 2013
 *      Author: ckow
 */

#include "BaldwinLomaxTurbulentViscosityStencil.h"

BaldwinLomaxTurbulentViscosityStencil::BaldwinLomaxTurbulentViscosityStencil(
		Parameters & parameters) :
		ComputeLocalViscosityStencil(parameters), _kappa(parameters.blm.kappa), _delta99(
				parameters.blm.delta99), _A_plus(26.0), _alpha(0.0168), _Ccp(
				1.6), _Cklep(0.3), _Cwk(0.25), _k(0.4) {
	// TODO Auto-generated constructor stub

}

BaldwinLomaxTurbulentViscosityStencil::~BaldwinLomaxTurbulentViscosityStencil() {
	// TODO Auto-generated destructor stub
}

