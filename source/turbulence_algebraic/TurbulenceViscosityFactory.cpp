/*
 * TurbulenceViscosityFactory.cpp
 *
 *  Created on: Dec 13, 2013
 *      Author: ckow
 */

#include "TurbulenceViscosityFactory.h"
#include "stencils/TurbulencePrandtlViscosityStencil.h"

TurbulenceViscosityFactory::TurbulenceViscosityFactory(Parameters & parameters) :
		_parameters(parameters) {
	// TODO Auto-generated constructor stub

}

TurbulenceViscosityFactory::~TurbulenceViscosityFactory() {
	// TODO Auto-generated destructor stub
}

ComputeLocalViscosityStencil* TurbulenceViscosityFactory::getLocalViscosityStencil() {
  ComputeLocalViscosityStencil* buffer = NULL;

  if (_parameters.blm.modelType == "prandtl") {
    buffer = new TurbulencePrandtlViscosityStencil(_parameters);
  }

  if (buffer == NULL){
    handleError(1,"buffer==NULL!");
  }
  return buffer;
}
