/*
 * TurbulenceViscosityFactory.h
 *
 *  Created on: Dec 13, 2013
 *      Author: ckow
 */

#ifndef TURBULENCEVISCOSITYFACTORY_H_
#define TURBULENCEVISCOSITYFACTORY_H_

#include "../Parameters.h"
#include "stencils/ComputeLocalViscosityStencil.h"

class TurbulenceViscosityFactory {
private:
	Parameters &_parameters;
public:
	TurbulenceViscosityFactory(Parameters & parameters);
	virtual ~TurbulenceViscosityFactory();

	ComputeLocalViscosityStencil* getLocalViscosityStencil();
};

#endif /* TURBULENCEVISCOSITYFACTORY_H_ */
