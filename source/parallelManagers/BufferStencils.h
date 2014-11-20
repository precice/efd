#ifndef _BUFFER_STENCILS_H_
#define _BUFFER_STENCILS_H_

#include "../FlowField.h"
#include "../Stencil.h"

// This class should not be called global

class VelocityBufferFillStencil : public BoundaryStencil<FlowField> {
private:
  FLOAT*& _leftBufferOut;
  FLOAT*& _rightBufferOut;
  FLOAT*& _bottomBufferOut;
  FLOAT*& _topBufferOut;
  FLOAT*& _frontBufferOut;
  FLOAT*& _backBufferOut;

public:
  VelocityBufferFillStencil(const Parameters& parameters,
                            FLOAT*&           leftBufferOut,
                            FLOAT*&           rightBufferOut,
                            FLOAT*&           bottomBufferOut,
                            FLOAT*&           topBufferOut,
                            FLOAT*&           frontBufferOut,
                            FLOAT*&           backBufferOut);

  void
  applyLeftWall(FlowField& flowField, int i, int j);

  void
  applyRightWall(FlowField& flowField, int i, int j);

  void
  applyBottomWall(FlowField& flowField, int i, int j);

  void
  applyTopWall(FlowField& flowField, int i, int j);

  void
  applyLeftWall(FlowField& flowField, int i, int j, int k);

  void
  applyRightWall(FlowField& flowField, int i, int j, int k);

  void
  applyBottomWall(FlowField& flowField, int i, int j, int k);

  void
  applyTopWall(FlowField& flowField, int i, int j, int k);

  void
  applyFrontWall(FlowField& flowField, int i, int j, int k);

  void
  applyBackWall(FlowField& flowField, int i, int j, int k);
};

class VelocityBufferReadStencil : public BoundaryStencil<FlowField> {
private:
  FLOAT*& _leftBufferIn;
  FLOAT*& _rightBufferIn;
  FLOAT*& _bottomBufferIn;
  FLOAT*& _topBufferIn;
  FLOAT*& _frontBufferIn;
  FLOAT*& _backBufferIn;

public:
  VelocityBufferReadStencil(const Parameters& parameters,
                            FLOAT*&           leftBufferIn,
                            FLOAT*&           rightBufferIn,
                            FLOAT*&           bottomBufferIn,
                            FLOAT*&           topBufferIn,
                            FLOAT*&           frontBufferIn,
                            FLOAT*&           backBufferIn);

  void
  applyLeftWall(FlowField& flowField, int i, int j);

  void
  applyRightWall(FlowField& flowField, int i, int j);

  void
  applyBottomWall(FlowField& flowField, int i, int j);

  void
  applyTopWall(FlowField& flowField, int i, int j);

  void
  applyLeftWall(FlowField& flowField, int i, int j, int k);

  void
  applyRightWall(FlowField& flowField, int i, int j, int k);

  void
  applyBottomWall(FlowField& flowField, int i, int j, int k);

  void
  applyTopWall(FlowField& flowField, int i, int j, int k);

  void
  applyFrontWall(FlowField& flowField, int i, int j, int k);

  void
  applyBackWall(FlowField& flowField, int i, int j, int k);
};

class PressureBufferFillStencil : public BoundaryStencil<FlowField> {
private:
  FLOAT*& _leftBufferOut;
  FLOAT*& _rightBufferOut;
  FLOAT*& _bottomBufferOut;
  FLOAT*& _topBufferOut;
  FLOAT*& _frontBufferOut;
  FLOAT*& _backBufferOut;

public:
  PressureBufferFillStencil(const Parameters& parameters,
                            FLOAT*&           leftBufferOut,
                            FLOAT*&           rightBufferOut,
                            FLOAT*&           bottomBufferOut,
                            FLOAT*&           topBufferOut,
                            FLOAT*&           frontBufferOut,
                            FLOAT*&           backBufferOut);

  void
  applyLeftWall(FlowField& flowField, int i, int j);

  void
  applyRightWall(FlowField& flowField, int i, int j);

  void
  applyBottomWall(FlowField& flowField, int i, int j);

  void
  applyTopWall(FlowField& flowField, int i, int j);

  void
  applyLeftWall(FlowField& flowField, int i, int j, int k);

  void
  applyRightWall(FlowField& flowField, int i, int j, int k);

  void
  applyBottomWall(FlowField& flowField, int i, int j, int k);

  void
  applyTopWall(FlowField& flowField, int i, int j, int k);

  void
  applyFrontWall(FlowField& flowField, int i, int j, int k);

  void
  applyBackWall(FlowField& flowField, int i, int j, int k);
};

class PressureBufferReadStencil : public BoundaryStencil<FlowField> {
private:
  FLOAT*& _leftBufferIn;
  FLOAT*& _rightBufferIn;
  FLOAT*& _bottomBufferIn;
  FLOAT*& _topBufferIn;
  FLOAT*& _frontBufferIn;
  FLOAT*& _backBufferIn;

public:
  PressureBufferReadStencil(const Parameters& parameters,
                            FLOAT*&           leftBufferIn,
                            FLOAT*&           rightBufferIn,
                            FLOAT*&           bottomBufferIn,
                            FLOAT*&           topBufferIn,
                            FLOAT*&           frontBufferIn,
                            FLOAT*&           backBufferIn);

  // Right now, they send an array bigger than necessary, corresponding to the
  // whole face of
  // the flow field, and not only those in the domain. It was easier to program
  // that way.

  void
  applyLeftWall(FlowField& flowField, int i, int j);

  void
  applyRightWall(FlowField& flowField, int i, int j);

  void
  applyBottomWall(FlowField& flowField, int i, int j);

  void
  applyTopWall(FlowField& flowField, int i, int j);

  void
  applyLeftWall(FlowField& flowField, int i, int j, int k);

  void
  applyRightWall(FlowField& flowField, int i, int j, int k);

  void
  applyBottomWall(FlowField& flowField, int i, int j, int k);

  void
  applyTopWall(FlowField& flowField, int i, int j, int k);

  void
  applyFrontWall(FlowField& flowField, int i, int j, int k);

  void
  applyBackWall(FlowField& flowField, int i, int j, int k);
};

#endif
