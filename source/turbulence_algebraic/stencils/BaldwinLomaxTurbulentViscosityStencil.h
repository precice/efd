#ifndef BALDWINLOMAXTURBULENTVISCOSITYSTENCIL_H_
#define BALDWINLOMAXTURBULENTVISCOSITYSTENCIL_H_

#include "ComputeLocalViscosityStencil.h"

class BaldwinLomaxTurbulentViscosityStencil : public
                                              ComputeLocalViscosityStencil {
private:
  FLOAT _kappa;
  FLOAT _delta99;
  FLOAT _A_plus;
  FLOAT _alpha;
  FLOAT _Ccp;
  FLOAT _Cklep;
  FLOAT _Cwk;
  FLOAT _k;

public:
  BaldwinLomaxTurbulentViscosityStencil(Parameters& parameters);

  virtual
  ~BaldwinLomaxTurbulentViscosityStencil();
};

#endif /* BALDWINLOMAXTURBULENTVISCOSITYSTENCIL_H_ */
