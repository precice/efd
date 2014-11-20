#ifndef _COMPUTELOCALVISCOSITYSTENCIL_H_
#define _COMPUTELOCALVISCOSITYSTENCIL_H_

#include "Stencil.h"
#include "TurbulentFlowField.h"

/** computes the local viscosity value.
 *  @author Philipp Neumann
 *  TODO DMITRII: Implement your things here
 */
class ComputeLocalViscosityStencil : public FieldStencil<TurbulentFlowField> {
protected:
  FLOAT _localVelocity[27 * 3];
  FLOAT _localMeshsize[27 * 3];

public:
  ComputeLocalViscosityStencil(const Parameters& parameters);

  virtual
  ~ComputeLocalViscosityStencil();
};

#endif // _COMPUTELOCALVISCOSITYSTENCIL_H_
