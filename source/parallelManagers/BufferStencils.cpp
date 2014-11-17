#include "BufferStencils.h"
VelocityBufferFillStencil::
VelocityBufferFillStencil(const Parameters& parameters,
                          FLOAT*&           leftBufferOut,
                          FLOAT*&           rightBufferOut,
                          FLOAT*&           bottomBufferOut,
                          FLOAT*&           topBufferOut,
                          FLOAT*&           frontBufferOut,
                          FLOAT*&           backBufferOut)
  : BoundaryStencil<FlowField>(parameters),
    _leftBufferOut(leftBufferOut),
    _rightBufferOut(rightBufferOut),
    _bottomBufferOut(bottomBufferOut),
    _topBufferOut(topBufferOut),
    _frontBufferOut(frontBufferOut),
    _backBufferOut(backBufferOut)
{}

void
VelocityBufferFillStencil::
applyLeftWall(FlowField& flowField, int i, int j) {
  _leftBufferOut[2 * j]     = flowField.getVelocity().getVector(2, j)[0];
  _leftBufferOut[2 * j + 1] = flowField.getVelocity().getVector(2, j)[1];
}

void
VelocityBufferFillStencil::
applyRightWall(FlowField& flowField, int i, int j) {
  _rightBufferOut[2 * j] = flowField.getVelocity().getVector(
    flowField.getNx(), j)[0];
  _rightBufferOut[2 * j + 1] = flowField.getVelocity().getVector(
    flowField.getNx() + 1, j)[1];
}

void
VelocityBufferFillStencil::
applyBottomWall(FlowField& flowField, int i, int j) {
  _bottomBufferOut[2 * i]     = flowField.getVelocity().getVector(i, 2)[0];
  _bottomBufferOut[2 * i + 1] = flowField.getVelocity().getVector(i, 2)[1];
}

void
VelocityBufferFillStencil::
applyTopWall(FlowField& flowField, int i, int j) {
  _topBufferOut[2 * i] = flowField.getVelocity().getVector(i,
                                                           flowField.getNy()
                                                           + 1)[0];
  _topBufferOut[2 * i + 1] = flowField.getVelocity().getVector(i,
                                                               flowField.getNy())
                             [1];
}

void
VelocityBufferFillStencil::
applyLeftWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (j + (flowField.getNy() + 3) * k);
  _leftBufferOut[index]     = flowField.getVelocity().getVector(2, j, k)[0];
  _leftBufferOut[index + 1] = flowField.getVelocity().getVector(2, j, k)[1];
  _leftBufferOut[index + 2] = flowField.getVelocity().getVector(2, j, k)[2];
}

void
VelocityBufferFillStencil::
applyRightWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (j + (flowField.getNy() + 3) * k);
  _rightBufferOut[index] = flowField.getVelocity().getVector(
    flowField.getNx(),   j, k)[0];
  _rightBufferOut[index + 1] = flowField.getVelocity().getVector(
    flowField.getNx() + 1, j, k)[1];
  _rightBufferOut[index + 2] = flowField.getVelocity().getVector(
    flowField.getNx() + 1, j, k)[2];
}

void
VelocityBufferFillStencil::
applyBottomWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (i + (flowField.getNx() + 3) * k);
  _bottomBufferOut[index]     = flowField.getVelocity().getVector(i, 2, k)[0];
  _bottomBufferOut[index + 1] = flowField.getVelocity().getVector(i, 2, k)[1];
  _bottomBufferOut[index + 2] = flowField.getVelocity().getVector(i, 2, k)[2];
}

void
VelocityBufferFillStencil::
applyTopWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (i + (flowField.getNx() + 3) * k);
  _topBufferOut[index] = flowField.getVelocity().getVector(i,
                                                           flowField.getNy()
                                                           + 1, k)[0];
  _topBufferOut[index + 1] = flowField.getVelocity().getVector(i,
                                                               flowField.getNy(),
                                                               k)[1];
  _topBufferOut[index + 2] = flowField.getVelocity().getVector(i,
                                                               flowField.getNy()
                                                               + 1, k)[2];
}

void
VelocityBufferFillStencil::
applyFrontWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (i + (flowField.getNx() + 3) * j);
  _frontBufferOut[index]     = flowField.getVelocity().getVector(i, j, 2)[0];
  _frontBufferOut[index + 1] = flowField.getVelocity().getVector(i, j, 2)[1];
  _frontBufferOut[index + 2] = flowField.getVelocity().getVector(i, j, 2)[2];
}

void
VelocityBufferFillStencil::
applyBackWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (i + (flowField.getNx() + 3) * j);
  _backBufferOut[index] = flowField.getVelocity().getVector(i, j,
                                                            flowField.getNz()
                                                            + 1)[0];
  _backBufferOut[index + 1] = flowField.getVelocity().getVector(i, j,
                                                                flowField.getNz()
                                                                + 1)[1];
  _backBufferOut[index + 2] = flowField.getVelocity().getVector(i, j,
                                                                flowField.getNz())
                              [2];
}

VelocityBufferReadStencil::
VelocityBufferReadStencil(const Parameters& parameters,
                          FLOAT*&           leftBufferIn,
                          FLOAT*&           rightBufferIn,
                          FLOAT*&           bottomBufferIn,
                          FLOAT*&           topBufferIn,
                          FLOAT*&           frontBufferIn,
                          FLOAT*&           backBufferIn)
  : BoundaryStencil<FlowField>(parameters),
    _leftBufferIn(leftBufferIn),
    _rightBufferIn(rightBufferIn),
    _bottomBufferIn(bottomBufferIn),
    _topBufferIn(topBufferIn),
    _frontBufferIn(frontBufferIn),
    _backBufferIn(backBufferIn)
{}

void
VelocityBufferReadStencil::
applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(0, j)[0] = _leftBufferIn[2 * j];
  flowField.getVelocity().getVector(1, j)[1] = _leftBufferIn[2 * j + 1];
}

void
VelocityBufferReadStencil::
applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(flowField.getNx() + 2, j)[0] =
    _rightBufferIn[2 * j];
  flowField.getVelocity().getVector(flowField.getNx() + 2, j)[1] =
    _rightBufferIn[2 * j + 1];
}

void
VelocityBufferReadStencil::
applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, 1)[0] = _bottomBufferIn[2 * i];
  flowField.getVelocity().getVector(i, 0)[1] = _bottomBufferIn[2 * i + 1];
}

void
VelocityBufferReadStencil::
applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, flowField.getNy() + 2)[0] =
    _topBufferIn[2 * i];
  flowField.getVelocity().getVector(i, flowField.getNy() + 2)[1] =
    _topBufferIn[2 * i + 1];
}

void
VelocityBufferReadStencil::
applyLeftWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (j + (flowField.getNy() + 3) * k);
  flowField.getVelocity().getVector(0, j, k)[0] = _leftBufferIn[index];
  flowField.getVelocity().getVector(1, j, k)[1] = _leftBufferIn[index + 1];
  flowField.getVelocity().getVector(1, j, k)[2] = _leftBufferIn[index + 2];
}

void
VelocityBufferReadStencil::
applyRightWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (j + (flowField.getNy() + 3) * k);
  flowField.getVelocity().getVector(flowField.getNx() + 2, j, k)[0] =
    _rightBufferIn[index];
  flowField.getVelocity().getVector(flowField.getNx() + 2, j, k)[1] =
    _rightBufferIn[index + 1];
  flowField.getVelocity().getVector(flowField.getNx() + 2, j, k)[2] =
    _rightBufferIn[index + 2];
}

void
VelocityBufferReadStencil::
applyBottomWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (i + (flowField.getNx() + 3) * k);
  flowField.getVelocity().getVector(i, 1, k)[0] = _bottomBufferIn[index];
  flowField.getVelocity().getVector(i, 0, k)[1] = _bottomBufferIn[index + 1];
  flowField.getVelocity().getVector(i, 1, k)[2] = _bottomBufferIn[index + 2];
}

void
VelocityBufferReadStencil::
applyTopWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (i + (flowField.getNx() + 3) * k);
  flowField.getVelocity().getVector(i, flowField.getNy() + 2, k)[0] =
    _topBufferIn[index];
  flowField.getVelocity().getVector(i, flowField.getNy() + 2, k)[1] =
    _topBufferIn[index + 1];
  flowField.getVelocity().getVector(i, flowField.getNy() + 2, k)[2] =
    _topBufferIn[index + 2];
}

void
VelocityBufferReadStencil::
applyFrontWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (i + (flowField.getNx() + 3) * j);
  flowField.getVelocity().getVector(i, j, 1)[0] = _frontBufferIn[index];
  flowField.getVelocity().getVector(i, j, 1)[1] = _frontBufferIn[index + 1];
  flowField.getVelocity().getVector(i, j, 0)[2] = _frontBufferIn[index + 2];
}

void
VelocityBufferReadStencil::
applyBackWall(FlowField& flowField, int i, int j, int k) {
  const int index = 3 * (i + (flowField.getNx() + 3) * j);
  flowField.getVelocity().getVector(i, j, flowField.getNz() + 2)[0] =
    _backBufferIn[index];
  flowField.getVelocity().getVector(i, j, flowField.getNz() + 2)[1] =
    _backBufferIn[index + 1];
  flowField.getVelocity().getVector(i, j, flowField.getNz() + 2)[2] =
    _backBufferIn[index + 2];
}

PressureBufferFillStencil::
PressureBufferFillStencil(const Parameters& parameters,
                          FLOAT*&           leftBufferOut,
                          FLOAT*&           rightBufferOut,
                          FLOAT*&           bottomBufferOut,
                          FLOAT*&           topBufferOut,
                          FLOAT*&           frontBufferOut,
                          FLOAT*&           backBufferOut)
  : BoundaryStencil<FlowField>(parameters),
    _leftBufferOut(leftBufferOut),
    _rightBufferOut(rightBufferOut),
    _bottomBufferOut(bottomBufferOut),
    _topBufferOut(topBufferOut),
    _frontBufferOut(frontBufferOut),
    _backBufferOut(backBufferOut)
{}

void
PressureBufferFillStencil::
applyLeftWall(FlowField& flowField, int i, int j) {
  _leftBufferOut[j] = flowField.getPressure().getScalar(2, j);
}

void
PressureBufferFillStencil::
applyRightWall(FlowField& flowField, int i, int j) {
  _rightBufferOut[j] = flowField.getPressure().getScalar(flowField.getNx() + 1,
                                                         j);
}

void
PressureBufferFillStencil::
applyBottomWall(FlowField& flowField, int i, int j) {
  _bottomBufferOut[i] = flowField.getPressure().getScalar(i, 2);
}

void
PressureBufferFillStencil::
applyTopWall(FlowField& flowField, int i, int j) {
  _topBufferOut[i] = flowField.getPressure().getScalar(i, flowField.getNy() +
                                                       1);
}

void
PressureBufferFillStencil::
applyLeftWall(FlowField& flowField, int i, int j, int k) {
  const int index = j + (flowField.getNy() + 3) * k;
  _leftBufferOut[index] = flowField.getPressure().getScalar(2, j, k);
}

void
PressureBufferFillStencil::
applyRightWall(FlowField& flowField, int i, int j, int k) {
  const int index = j + (flowField.getNy() + 3) * k;
  _rightBufferOut[index] = flowField.getPressure().getScalar(flowField.getNx() +
                                                             1, j, k);
}

void
PressureBufferFillStencil::
applyBottomWall(FlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * k;
  _bottomBufferOut[index] = flowField.getPressure().getScalar(i, 2, k);
}

void
PressureBufferFillStencil::
applyTopWall(FlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * k;
  _topBufferOut[index] = flowField.getPressure().getScalar(i,
                                                           flowField.getNy() +
                                                           1, k);
}

void
PressureBufferFillStencil::
applyFrontWall(FlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * j;
  _frontBufferOut[index] = flowField.getPressure().getScalar(i, j, 2);
}

void
PressureBufferFillStencil::
applyBackWall(FlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * j;
  _backBufferOut[index] = flowField.getPressure().getScalar(i, j,
                                                            flowField.getNz() +
                                                            1);
}

PressureBufferReadStencil::
PressureBufferReadStencil(const Parameters& parameters,
                          FLOAT*&           leftBufferIn,
                          FLOAT*&           rightBufferIn,
                          FLOAT*&           bottomBufferIn,
                          FLOAT*&           topBufferIn,
                          FLOAT*&           frontBufferIn,
                          FLOAT*&           backBufferIn)
  : BoundaryStencil<FlowField>(parameters),
    _leftBufferIn(leftBufferIn),
    _rightBufferIn(rightBufferIn),
    _bottomBufferIn(bottomBufferIn),
    _topBufferIn(topBufferIn),
    _frontBufferIn(frontBufferIn),
    _backBufferIn(backBufferIn)
{}

void
PressureBufferReadStencil::
applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(1, j) = _leftBufferIn[j];
}

void
PressureBufferReadStencil::
applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(flowField.getNx() + 2, j) =
    _rightBufferIn[j];
}

void
PressureBufferReadStencil::
applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, 1) = _bottomBufferIn[i];
}

void
PressureBufferReadStencil::
applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, flowField.getNy() + 2) = _topBufferIn[i];
}

void
PressureBufferReadStencil::
applyLeftWall(FlowField& flowField, int i, int j, int k) {
  const int index = j + (flowField.getNy() + 3) * k;
  flowField.getPressure().getScalar(1, j, k) = _leftBufferIn[index];
}

void
PressureBufferReadStencil::
applyRightWall(FlowField& flowField, int i, int j, int k) {
  const int index = j + (flowField.getNy() + 3) * k;
  flowField.getPressure().getScalar(flowField.getNx() + 2, j, k) =
    _rightBufferIn[index];
}

void
PressureBufferReadStencil::
applyBottomWall(FlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * k;
  flowField.getPressure().getScalar(i, 1, k) = _bottomBufferIn[index];
}

void
PressureBufferReadStencil::
applyTopWall(FlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * k;
  flowField.getPressure().getScalar(i, flowField.getNy() + 2, k) =
    _topBufferIn[index];
}

void
PressureBufferReadStencil::
applyFrontWall(FlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * j;
  flowField.getPressure().getScalar(i, j, 1) = _frontBufferIn[index];
}

void
PressureBufferReadStencil::
applyBackWall(FlowField& flowField, int i, int j, int k) {
  const int index = i + (flowField.getNx() + 3) * j;
  flowField.getPressure().getScalar(i, j, flowField.getNz() + 2) =
    _backBufferIn[index];
}
