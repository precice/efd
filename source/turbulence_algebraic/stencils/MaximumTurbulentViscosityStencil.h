#ifndef _MAXIMUMTURBULENTVISCOSITYSTENCIL_H_
#define _MAXIMUMTURBULENTVISCOSITYSTENCIL_H_

#include "../Stencil.h"
#include "../TurbulentFlowField.h"

/** compute the maximum turbulent viscosity, required for time step/stability
 * restrictions. For stretched meshes, we currently compute the max. turbulent
 * viscosity and evaluate the stability criterion together with the min.
 * meshsize throughout the whole domain.
 * TODO: improve this and evaluate maximum value for local meshsize+local
 * viscosity.
 *  @author Philipp Neumann
 */
class MaximumTurbulentViscosityStencil : public
                                         FieldStencil<TurbulentFlowField> {
private:
  FLOAT _maximumTurbulentViscosity;

public:
  MaximumTurbulentViscosityStencil(const Parameters& parameters)
    : FieldStencil<TurbulentFlowField>(parameters), _maximumTurbulentViscosity(
        0.0) {}
  virtual
  ~MaximumTurbulentViscosityStencil() {}

  void
  reset() {
    _maximumTurbulentViscosity = 0.0;
  }
  FLOAT
  getMaximumTurbulentViscosity() const { return _maximumTurbulentViscosity; }

  void
  apply(TurbulentFlowField& flowField, int i, int j) {
    if (flowField.getViscosityField().getScalar(i, j) >
        _maximumTurbulentViscosity) {
      _maximumTurbulentViscosity = flowField.getViscosityField().getScalar(i,
                                                                           j);
    }
  }
  void
  apply(TurbulentFlowField& flowField, int i, int j, int k) {
    if (flowField.getViscosityField().getScalar(i, j, k) >
        _maximumTurbulentViscosity) {
      _maximumTurbulentViscosity = flowField.getViscosityField().getScalar(i, j,
                                                                           k);
    }
  }
};

#endif // _MAXIMUMTURBULENTVISCOSITYSTENCIL_H_
