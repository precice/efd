#include "../stencils/TurbulenceStencilFunctions.h"
#include "ComputeLocalViscosityStencil.h"
//

ComputeLocalViscosityStencil::
ComputeLocalViscosityStencil(
  const Parameters& parameters)
  : FieldStencil<TurbulentFlowField>(parameters) {}

ComputeLocalViscosityStencil::~ComputeLocalViscosityStencil() {}
