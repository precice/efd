#include "ComputeLocalViscosityStencil.h"
#include "../stencils/TurbulenceStencilFunctions.h"

ComputeLocalViscosityStencil::ComputeLocalViscosityStencil(
		const Parameters &parameters) :
		FieldStencil<TurbulentFlowField>(parameters){}

ComputeLocalViscosityStencil::~ComputeLocalViscosityStencil() {}

