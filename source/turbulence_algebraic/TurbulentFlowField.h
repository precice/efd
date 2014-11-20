#ifndef _TURBULENTFLOWFIELD_H_
#define _TURBULENTFLOWFIELD_H_

#include "../FlowField.h"

/** data structure for flow field in case of turbulent simulation.
 *  @author Philipp Neumann
 */
class TurbulentFlowField : public FlowField {
private:
  /** scalar field for viscosity values */
  ScalarField _viscosity;
  ScalarField _distanceToWall;

public:
  /** initialise the default flow field AND the new data field for viscosity */
  TurbulentFlowField(const Parameters& parameters);

  ScalarField&
  getViscosityField();

  ScalarField&
  getDistanceToWallField();
};

#endif // _TURBULENTFLOWFIELD_H_
