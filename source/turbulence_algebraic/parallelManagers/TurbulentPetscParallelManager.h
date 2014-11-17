#ifndef _TURBULENTPETSCPARALLELMANAGER_H_
#define _TURBULENTPETSCPARALLELMANAGER_H_

#include "parallelManagers/PetscParallelManager.h"
#include "parallelManagers/FillReadViscosityBufferStencil.h"


/** extends the PetscParallelManager by communication of the viscosity values across parallel boundaries.
 *  The whole exchange mechanism works completely the same as for the pressure, but reading/writing viscosity values from
 *  the new data structure TurbulentFlowField instead.
 *  @author Philipp Neumann
 */
class TurbulentPetscParallelManager: public PetscParallelManager {
  public:
    TurbulentPetscParallelManager (TurbulentFlowField & flowField, const Parameters & parameters);
    virtual ~TurbulentPetscParallelManager();

    void communicateViscosity();

  private:
    ViscosityBufferFillStencil _fillViscosityStencil;
    ViscosityBufferReadStencil _readViscosityStencil;
    ParallelBoundaryIterator<TurbulentFlowField> _fillViscosityIterator;
    ParallelBoundaryIterator<TurbulentFlowField> _readViscosityIterator;
};

#endif // _TURBULENTPETSCPARALLELMANAGER_H_
