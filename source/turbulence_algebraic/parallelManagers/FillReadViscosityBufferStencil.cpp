#include "parallelManagers/FillReadViscosityBufferStencil.h"
//

ViscosityBufferFillStencil::
ViscosityBufferFillStencil(const Parameters& parameters,
                           FLOAT*&           leftBufferOut,
                           FLOAT*&           rightBufferOut,
                           FLOAT*&           bottomBufferOut,
                           FLOAT*&           topBufferOut,
                           FLOAT*&           frontBufferOut,
                           FLOAT*&           backBufferOut)
  : BoundaryStencil<TurbulentFlowField>(parameters),
    _leftBufferOut(leftBufferOut),
    _rightBufferOut(rightBufferOut),
    _bottomBufferOut(bottomBufferOut),
    _topBufferOut(topBufferOut),
    _frontBufferOut(frontBufferOut),
    _backBufferOut(backBufferOut)
{}

void
ViscosityBufferFillStencil::
applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  _leftBufferOut[j] = flowField.getViscosityField().getScalar(2, j);
}

void
ViscosityBufferFillStencil::
applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  _rightBufferOut[j] = flowField.getViscosityField().getScalar(
    flowField.getNx() + 1, j);
}

void
ViscosityBufferFillStencil::
applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  _bottomBufferOut[i] = flowField.getViscosityField().getScalar(i, 2);
}

void
ViscosityBufferFillStencil::
applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  _topBufferOut[i] = flowField.getViscosityField().getScalar(i,
                                                             flowField.getNy() +
                                                             1);
}

void
ViscosityBufferFillStencil::
applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = j + (flowField.getNy() + 3) * k;
  _leftBufferOut[index] = flowField.getViscosityField().getScalar(2, j, k);
}

void
ViscosityBufferFillStencil::
applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = j + (flowField.getNy() + 3) * k;
  _rightBufferOut[index] = flowField.getViscosityField().getScalar(
    flowField.getNx() + 1, j, k);
}

void
ViscosityBufferFillStencil::
applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * k;
  _bottomBufferOut[index] = flowField.getViscosityField().getScalar(i, 2, k);
}

void
ViscosityBufferFillStencil::
applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * k;
  _topBufferOut[index] = flowField.getViscosityField().getScalar(i,
                                                                 flowField.getNy()
                                                                 + 1, k);
}

void
ViscosityBufferFillStencil::
applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * j;
  _frontBufferOut[index] = flowField.getViscosityField().getScalar(i, j, 2);
}

void
ViscosityBufferFillStencil::
applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * j;
  _backBufferOut[index] = flowField.getViscosityField().getScalar(i, j,
                                                                  flowField.
                                                                  getNz() + 1);
}

ViscosityBufferReadStencil::
ViscosityBufferReadStencil(const Parameters& parameters,
                           FLOAT*&           leftBufferIn,
                           FLOAT*&           rightBufferIn,
                           FLOAT*&           bottomBufferIn,
                           FLOAT*&           topBufferIn,
                           FLOAT*&           frontBufferIn,
                           FLOAT*&           backBufferIn)
  : BoundaryStencil<TurbulentFlowField>(parameters),
    _leftBufferIn(leftBufferIn),
    _rightBufferIn(rightBufferIn),
    _bottomBufferIn(bottomBufferIn),
    _topBufferIn(topBufferIn),
    _frontBufferIn(frontBufferIn),
    _backBufferIn(backBufferIn)
{}

void
ViscosityBufferReadStencil::
applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosityField().getScalar(1, j) = _leftBufferIn[j];
}

void
ViscosityBufferReadStencil::
applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosityField().getScalar(flowField.getNx() + 2, j) =
    _rightBufferIn[j];
}

void
ViscosityBufferReadStencil::
applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosityField().getScalar(i, 1) = _bottomBufferIn[i];
}

void
ViscosityBufferReadStencil::
applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosityField().getScalar(i, flowField.getNy() + 2) =
    _topBufferIn[i];
}

void
ViscosityBufferReadStencil::
applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = j + (flowField.getNy() + 3) * k;
  flowField.getViscosityField().getScalar(1, j, k) = _leftBufferIn[index];
}

void
ViscosityBufferReadStencil::
applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = j + (flowField.getNy() + 3) * k;
  flowField.getViscosityField().getScalar(flowField.getNx() + 2, j, k) =
    _rightBufferIn[index];
}

void
ViscosityBufferReadStencil::
applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * k;
  flowField.getViscosityField().getScalar(i, 1, k) = _bottomBufferIn[index];
}

void
ViscosityBufferReadStencil::
applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * k;
  flowField.getViscosityField().getScalar(i, flowField.getNy() + 2, k) =
    _topBufferIn[index];
}

void
ViscosityBufferReadStencil::
applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * j;
  flowField.getViscosityField().getScalar(i, j, 1) = _frontBufferIn[index];
}

void
ViscosityBufferReadStencil::
applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * j;
  flowField.getViscosityField().getScalar(i, j, flowField.getNz() + 2) =
    _backBufferIn[index];
}
