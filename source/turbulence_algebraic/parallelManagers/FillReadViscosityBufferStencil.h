#include "Stencil.h"
#include "TurbulentFlowField.h"

class ViscosityBufferFillStencil : public BoundaryStencil<TurbulentFlowField> {
private:
  FLOAT*& _leftBufferOut;
  FLOAT*& _rightBufferOut;
  FLOAT*& _bottomBufferOut;
  FLOAT*& _topBufferOut;
  FLOAT*& _frontBufferOut;
  FLOAT*& _backBufferOut;

public:
  ViscosityBufferFillStencil(const Parameters& parameters,
                             FLOAT*&           leftBufferOut,
                             FLOAT*&           rightBufferOut,
                             FLOAT*&           bottomBufferOut,
                             FLOAT*&           topBufferOut,
                             FLOAT*&           frontBufferOut,
                             FLOAT*&           backBufferOut);

  void
  applyLeftWall(TurbulentFlowField& flowField, int i, int j);
  void
  applyRightWall(TurbulentFlowField& flowField, int i, int j);
  void
  applyBottomWall(TurbulentFlowField& flowField, int i, int j);
  void
  applyTopWall(TurbulentFlowField& flowField, int i, int j);

  void
  applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyRightWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyTopWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyBackWall(TurbulentFlowField& flowField, int i, int j, int k);
};

class ViscosityBufferReadStencil : public BoundaryStencil<TurbulentFlowField> {
private:
  FLOAT*& _leftBufferIn;
  FLOAT*& _rightBufferIn;
  FLOAT*& _bottomBufferIn;
  FLOAT*& _topBufferIn;
  FLOAT*& _frontBufferIn;
  FLOAT*& _backBufferIn;

public:
  ViscosityBufferReadStencil(const Parameters& parameters,
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
  applyLeftWall(TurbulentFlowField& flowField, int i, int j);
  void
  applyRightWall(TurbulentFlowField& flowField, int i, int j);
  void
  applyBottomWall(TurbulentFlowField& flowField, int i, int j);
  void
  applyTopWall(TurbulentFlowField& flowField, int i, int j);

  void
  applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyRightWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyTopWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k);
  void
  applyBackWall(TurbulentFlowField& flowField, int i, int j, int k);
};
