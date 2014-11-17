#include "TurbulentFlowField.h"

TurbulentFlowField::TurbulentFlowField(const Parameters & parameters) :
		FlowField(parameters), _viscosity(
				parameters.geometry.dim == 2 ?
						ScalarField(parameters.parallel.localSize[0] + 3,
								parameters.parallel.localSize[1] + 3) :
						ScalarField(parameters.parallel.localSize[0] + 3,
								parameters.parallel.localSize[1] + 3,
								parameters.parallel.localSize[2] + 3)), _distanceToWall(
				parameters.geometry.dim == 2 ?
						ScalarField(parameters.parallel.localSize[0] + 3,
								parameters.parallel.localSize[1] + 3) :
						ScalarField(parameters.parallel.localSize[0] + 3,
								parameters.parallel.localSize[1] + 3,
								parameters.parallel.localSize[2] + 3)) {
}

ScalarField& TurbulentFlowField::getViscosityField() {
	return _viscosity;
}

ScalarField& TurbulentFlowField::getDistanceToWallField(){
	return _distanceToWall;
}
