/*
 * TurbFGHStencil.h
 *
 *  Created on: Dec 13, 2013
 *      Author: ckow
 */

#ifndef TURBFGHSTENCIL_H_
#define TURBFGHSTENCIL_H_

#include "../Stencil.h"
#include "../TurbulentFlowField.h"

class TurbFGHStencil: public FieldStencil<TurbulentFlowField> {
private:
	FLOAT _localVelocity[27*3];
        FLOAT _localMeshsize[27*3];
	FLOAT _localViscosity[27];
public:
	TurbFGHStencil(const Parameters & parameters );

	void apply(TurbulentFlowField & flowField,int i,int j);
	void apply(TurbulentFlowField & flowField,int i,int j,int k);
};

#endif /* TURBFGHSTENCIL_H_ */
