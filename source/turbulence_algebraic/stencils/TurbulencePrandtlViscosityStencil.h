/*
 * TurbulencePrandtlViscosityStencil.h
 *
 *  Created on: Dec 13, 2013
 *      Author: ckow
 */

#ifndef TURBULENCEPRANDTLVISCOSITYSTENCIL_H_
#define TURBULENCEPRANDTLVISCOSITYSTENCIL_H_

#include "ComputeLocalViscosityStencil.h"
#include "TurbulentFlowField.h"

class TurbulencePrandtlViscosityStencil: public ComputeLocalViscosityStencil {
private:
	FLOAT _kappa;
	FLOAT _delta99;
public:
	TurbulencePrandtlViscosityStencil(Parameters &parameters);
	virtual ~TurbulencePrandtlViscosityStencil();
    void apply ( TurbulentFlowField & flowField, int i, int j);

    /** Performs the operation in 3D in a given position
     * @param turbulentflowField Flow field data for a turbulence simulation
     * @param parameters Parameters of the problem
     * @param i Position in the x direction
     * @param j Position in the y direction
     * @param k Position in the z direction
     */
    void apply ( TurbulentFlowField & flowField, int i, int j, int k);
};

#endif /* TURBULENCEPRANDTLVISCOSITYSTENCIL_H_ */
