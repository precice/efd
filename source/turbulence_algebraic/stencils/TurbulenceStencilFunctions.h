/*
 * TurbulenceStencilFunctions.h
 *
 *  Created on: Dec 12, 2013
 *      Author: kowitz
 */

#ifndef TURBULENCESTENCILFUNCTIONS_H_
#define TURBULENCESTENCILFUNCTIONS_H_

#include "../../stencils/StencilFunctions.h"

inline int mapds(int i,int j, int k) {
    return 13 + 9 * k + 3 * j + i;
}


// needed for F, see Griebel book, Eq. (10,28), modified for non-uniform meshes
inline FLOAT dx_du_dx(const FLOAT * const lv, const FLOAT * const lm, const FLOAT * const turbVisc, const Parameters & parameters) {
/*  const FLOAT tmp1= (  (turbVisc[mapds(1,0,0)]+1. / parameters.flow.Re)*(lv[mapd(1,0,0,0)] - lv[mapd(0,0,0,0)])
          - (turbVisc[mapds(0,0,0)]+1. / parameters.flow.Re)*(lv[mapd(0,0,0,0)] - lv[mapd(-1,0,0,0)])
         )/(lm[mapd(0,0,0,0)]*lm[mapd(0,0,0,0)]);
*/

  // first, we evaluate the expression visc*du/dx, using a central difference, in the middle of the right and current cell (current cell
  // corresponds to leftDifference)
  const FLOAT rightDifference = (lv[mapd(1,0,0,0)] - lv[mapd(0,0,0,0)])* (turbVisc[mapds(1,0,0)]+ 1.0/parameters.flow.Re)/lm[mapd(1,0,0,0)];
  const FLOAT leftDifference  = (lv[mapd(0,0,0,0)] - lv[mapd(-1,0,0,0)])*(turbVisc[mapds(0,0,0)]+ 1.0/parameters.flow.Re)/lm[mapd(0,0,0,0)];
  // next we compute a central difference from these two expressions to compute the second-order-like derivate. For uniform meshes,
  // this yield second-order accuracy.
  const FLOAT tmp2 = (rightDifference-leftDifference)/(0.5*(lm[mapd(0,0,0,0)]+lm[mapd(1,0,0,0)]));

//  if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dx_du_dx"); }
  return tmp2;
}


// needed for G, see Griebel book, Eq. (10,28). For details on implementation, see dx_du_dx
inline FLOAT dy_dv_dy(const FLOAT * const lv, const FLOAT * const lm, const FLOAT * const turbVisc, const Parameters & parameters) {
/*  const FLOAT tmp1= (  (turbVisc[mapds(0,1,0)]+1. / parameters.flow.Re)*(lv[mapd(0,1,0,1)] - lv[mapd(0,0,0,1)])
          - (turbVisc[mapds(0,0,0)]+1. / parameters.flow.Re)*(lv[mapd(0,0,0,1)] - lv[mapd(0,-1,0,1)])
         )/(lm[mapd(0,0,0,1)]*lm[mapd(0,0,0,1)]);
*/
  const FLOAT rightDifference = (lv[mapd(0,1,0,1)] - lv[mapd(0,0,0,1)])* (turbVisc[mapds(0,1,0)]+ 1.0/parameters.flow.Re)/lm[mapd(0,1,0,1)];
  const FLOAT leftDifference  = (lv[mapd(0,0,0,1)] - lv[mapd(0,-1,0,1)])*(turbVisc[mapds(0,0,0)]+ 1.0/parameters.flow.Re)/lm[mapd(0,0,0,1)];
  const FLOAT tmp2 = (rightDifference-leftDifference)/(0.5*(lm[mapd(0,0,0,1)]+lm[mapd(0,1,0,1)]));

//  if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dy_dv_dy"); }
  return tmp2;
}

// needed for H, see Griebel book, Eq. (10,28). For details on implementation, see dx_du_dx
inline FLOAT dz_dw_dz(const FLOAT * const lv, const FLOAT * const lm, const FLOAT * const turbVisc, const Parameters & parameters) {
/**  const FLOAT tmp1= (  (turbVisc[mapds(0,0,1)]+1. / parameters.flow.Re)*(lv[mapd(0,0,1,2)] - lv[mapd(0,0,0,2)])
          - (turbVisc[mapds(0,0,0)]+1. / parameters.flow.Re)*(lv[mapd(0,0,0,2)] - lv[mapd(0,0,-1,2)])
         )/(lm[mapd(0,0,0,2)]*lm[mapd(0,0,0,2)]);
*/
  const FLOAT rightDifference = (lv[mapd(0,0,1,2)] - lv[mapd(0,0,0,2)])* (turbVisc[mapds(0,0,1)]+ 1.0/parameters.flow.Re)/lm[mapd(0,0,1,2)];
  const FLOAT leftDifference  = (lv[mapd(0,0,0,2)] - lv[mapd(0,0,-1,2)])*(turbVisc[mapds(0,0,0)]+ 1.0/parameters.flow.Re)/lm[mapd(0,0,0,2)];
  const FLOAT tmp2 = (rightDifference-leftDifference)/(0.5*(lm[mapd(0,0,0,2)]+lm[mapd(0,0,1,2)]));

//  if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dz_dw_dz"); }
  return tmp2;
}


// carries out bilinear interpolation assuming the values visc00 - visc11 located at the points (0,0) to (lx,ly). The evaluation is carried out at the point (x,y)
inline FLOAT bilinear(const FLOAT visc00, const FLOAT visc10, const FLOAT visc01, const FLOAT visc11, const FLOAT lx, const FLOAT ly, const FLOAT x, const FLOAT y){
  return ( visc00*(lx-x)*(ly-y) + visc10*x*(ly-y) + visc01*(lx-x)*y + visc11*x*y)/(lx*ly);
}


// needed for F, see Griebel book, Eq. (10,28)
inline FLOAT dy_du_dy_plus_dv_dx(const FLOAT * const lv, const FLOAT * const lm, const FLOAT * const turbVisc, const Parameters & parameters) {
/*  const FLOAT tmp1visc1 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 1,0,0)]+turbVisc[mapds(0, 1,0)]+turbVisc[mapds( 1, 1,0)])+1. / parameters.flow.Re;
  const FLOAT tmp1visc2 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 1,0,0)]+turbVisc[mapds(0,-1,0)]+turbVisc[mapds( 1,-1,0)])+1. / parameters.flow.Re;

  const FLOAT tmp1= (tmp1visc1* (  (lv[mapd(0,1,0,0)] - lv[mapd(0,0,0,0)] )/lm[mapd(0,0,0,1)]
                                 + (lv[mapd(1,0,0,1)] - lv[mapd(0,0,0,1)] )/lm[mapd(0,0,0,0)] )
                    -tmp1visc2* (  (lv[mapd(0,0,0,0)] - lv[mapd(0,-1,0,0)])/lm[mapd(0,0,0,1)]
                                 + (lv[mapd(1,-1,0,1)]- lv[mapd(0,-1,0,1)])/lm[mapd(0,0,0,0)] ) )/lm[mapd(0,0,0,1)];
*/


  const FLOAT Lx1 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd(1,0,0,0)]);
  const FLOAT Ly0 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd(0,-1,0,1)]);
  const FLOAT Ly1 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd(0, 1,0,1)]);

  // first, we interpolate the viscosity values at (for 2D) lower right corner and upper right corner of current cell.
  const FLOAT visc0 = bilinear(turbVisc[mapds( 0,-1, 0)], turbVisc[mapds( 1,-1, 0)], turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 1, 0, 0)],
                               Lx1,Ly0,0.5*lm[mapd( 0,-1, 0, 0)],0.5*lm[mapd( 0,-1, 0, 1)]) + 1.0/parameters.flow.Re;
  const FLOAT visc1 = bilinear(turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 1, 0, 0)], turbVisc[mapds( 0, 1, 0)], turbVisc[mapds( 1, 1, 0)],
                               Lx1,Ly1,0.5*lm[mapd( 0, 0, 0, 0)],0.5*lm[mapd( 0, 0, 0, 1)]) + 1.0/parameters.flow.Re;
//  if (fabs(visc0-tmp1visc2)>1.0e-12){handleError(1,"dy_du_dy_plus_dv_dx_visc0");}
//  if (fabs(visc1-tmp1visc1)>1.0e-12){handleError(1,"dy_du_dy_plus_dv_dx_visc1");}

  // next, we evaluate the gradients du_dy and dv_dx at the lower right and upper right corners
  const FLOAT du_dy0 = (lv[mapd(0,0,0,0)] - lv[mapd(0,-1,0,0)])/Ly0;
  const FLOAT du_dy1 = (lv[mapd(0,1,0,0)] - lv[mapd(0, 0,0,0)])/Ly1;

  const FLOAT dv_dx0 = (lv[mapd(1,-1,0,1)]- lv[mapd(0,-1,0,1)])/Lx1;
  const FLOAT dv_dx1 = (lv[mapd(1, 0,0,1)]- lv[mapd(0, 0,0,1)])/Lx1;

  // finally, we compute the gradient in y-direction using a central difference along right cell edge (for 2D)
  const FLOAT tmp2 = ( visc1*(du_dy1+dv_dx1) - visc0*(du_dy0+dv_dx0) )/lm[mapd(0,0,0,1)];
//  if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dy_du_dy_plus_dv_dx");}
  return tmp2;
}


// for details on implementation, see dy_du_dy_plus_dv_dx
inline FLOAT dz_du_dz_plus_dw_dx(const FLOAT * const lv, const FLOAT * const lm, const FLOAT * const turbVisc, const Parameters & parameters) {
/*  const FLOAT tmp1visc1 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 1,0,0)]+turbVisc[mapds(0, 0, 1)]+turbVisc[mapds( 1, 0, 1)])+1. / parameters.flow.Re;
  const FLOAT tmp1visc2 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 1,0,0)]+turbVisc[mapds(0, 0,-1)]+turbVisc[mapds( 1, 0,-1)])+1. / parameters.flow.Re;

  const FLOAT tmp1 =  (tmp1visc1* (  (lv[mapd(0,0,1,0)] - lv[mapd(0,0,0,0)] )/lm[mapd(0,0,0,2)]
                                + (lv[mapd(1,0,0,2)] - lv[mapd(0,0,0,2)] )/lm[mapd(0,0,0,0)] )
                      -tmp1visc2* (  (lv[mapd(0,0,0,0)] - lv[mapd(0,0,-1,0)])/lm[mapd(0,0,0,2)]
                                   + (lv[mapd(1,0,-1,2)]- lv[mapd(0,0,-1,2)])/lm[mapd(0,0,0,0)] ) )/lm[mapd(0,0,0,2)];
*/

  const FLOAT Lx1 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd(1,0,0,0)]);
  const FLOAT Lz0 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd(0,0,-1,2)]);
  const FLOAT Lz1 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd(0,0, 1,2)]);

  const FLOAT visc0 = bilinear(turbVisc[mapds( 0, 0,-1)], turbVisc[mapds( 1, 0,-1)], turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 1, 0, 0)],
                               Lx1,Lz0,0.5*lm[mapd( 0, 0,-1, 0)],0.5*lm[mapd( 0, 0,-1, 2)]) + 1.0/parameters.flow.Re;
  const FLOAT visc1 = bilinear(turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 1, 0, 0)], turbVisc[mapds( 0, 0, 1)], turbVisc[mapds( 1, 0, 1)],
                               Lx1,Lz1,0.5*lm[mapd( 0, 0, 0, 0)],0.5*lm[mapd( 0, 0, 0, 2)]) + 1.0/parameters.flow.Re;
//  if (fabs(visc0-tmp1visc2)>1.0e-12){handleError(1,"dz_du_dz_plus_dw_dx_visc0");}
//  if (fabs(visc1-tmp1visc1)>1.0e-12){handleError(1,"dz_du_dz_plus_dw_dx_visc1");}

  const FLOAT du_dz0 = (lv[mapd(0,0,0,0)] - lv[mapd(0,0,-1,0)])/Lz0;
  const FLOAT du_dz1 = (lv[mapd(0,0,1,0)] - lv[mapd(0,0, 0,0)])/Lz1;

  const FLOAT dw_dx0 = (lv[mapd(1,0,-1,2)]- lv[mapd(0,0,-1,2)])/Lx1;
  const FLOAT dw_dx1 = (lv[mapd(1,0, 0,2)]- lv[mapd(0,0, 0,2)])/Lx1;

  const FLOAT tmp2 = ( visc1*(du_dz1+dw_dx1) - visc0*(du_dz0+dw_dx0) )/lm[mapd(0,0,0,2)];
//  if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dz_du_dx_plus_dw_dx");}
  return tmp2;
}


// needed for G, see Griebel book, Eq. (10,28); for details on implementation, see dy_du_dy_plus_dv_dx
inline FLOAT dx_dv_dx_plus_du_dy(const FLOAT * const lv, const FLOAT * const lm, const FLOAT * const turbVisc, const Parameters & parameters) {
/*  const FLOAT tmp1visc1 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 1,0,0)]+turbVisc[mapds(0, 1,0)]+turbVisc[mapds( 1, 1,0)])+1. / parameters.flow.Re;
  const FLOAT tmp1visc2 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 0,1,0)]+turbVisc[mapds(-1,0,0)]+turbVisc[mapds(-1, 1,0)])+1. / parameters.flow.Re;

  const FLOAT tmp1= (tmp1visc1* (  (lv[mapd(1,0,0,1)] - lv[mapd(0,0,0,1)] )/lm[mapd(0,0,0,0)]
                                 + (lv[mapd(0,1,0,0)] - lv[mapd(0,0,0,0)] )/lm[mapd(0,0,0,1)] )
                    -tmp1visc2* (  (lv[mapd(0,0,0,1)] - lv[mapd(-1,0,0,1)])/lm[mapd(0,0,0,0)]
                                 + (lv[mapd(-1,1,0,0)]- lv[mapd(-1,0,0,0)])/lm[mapd(0,0,0,1)] ) )/lm[mapd(0,0,0,0)];
*/

  const FLOAT Ly1 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd(0 ,1,0,1)]);
  const FLOAT Lx0 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd(-1,0,0,0)]);
  const FLOAT Lx1 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd( 1,0,0,0)]);

  const FLOAT visc0 = bilinear(turbVisc[mapds(-1, 0, 0)], turbVisc[mapds( 0, 0, 0)], turbVisc[mapds(-1, 1, 0)], turbVisc[mapds( 0, 1, 0)],
                               Lx0,Ly1,0.5*lm[mapd(-1, 0, 0, 0)],0.5*lm[mapd(-1, 0, 0, 1)]) + 1.0/parameters.flow.Re;
  const FLOAT visc1 = bilinear(turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 1, 0, 0)], turbVisc[mapds( 0, 1, 0)], turbVisc[mapds( 1, 1, 0)],
                               Lx1,Ly1,0.5*lm[mapd( 0, 0, 0, 0)],0.5*lm[mapd( 0, 0, 0, 1)]) + 1.0/parameters.flow.Re;
//  if (fabs(visc0-tmp1visc2)>1.0e-12){handleError(1,"dx_dv_dx_plus_du_dy_visc0");}
//  if (fabs(visc1-tmp1visc1)>1.0e-12){handleError(1,"dx_dv_dx_plus_du_dy_visc1");}

  const FLOAT du_dy0 = (lv[mapd(-1,1,0,0)] - lv[mapd(-1,0,0,0)])/Ly1;
  const FLOAT du_dy1 = (lv[mapd( 0,1,0,0)] - lv[mapd( 0,0,0,0)])/Ly1;

  const FLOAT dv_dx0 = (lv[mapd(0, 0,0,1)]- lv[mapd(-1,0,0,1)])/Lx0;
  const FLOAT dv_dx1 = (lv[mapd(1, 0,0,1)]- lv[mapd( 0,0,0,1)])/Lx1;

  const FLOAT tmp2 = ( visc1*(du_dy1+dv_dx1) - visc0*(du_dy0+dv_dx0) )/lm[mapd(0,0,0,0)];
//  if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dx_dv_dx_plus_du_dy");}
  return tmp2;
}


// for details on implementation, see dy_du_dy_plus_dv_dx
inline FLOAT dz_dv_dz_plus_dw_dy(const FLOAT * const lv, const FLOAT * const lm, const FLOAT * const turbVisc, const Parameters & parameters) {
/*  const FLOAT tmp1visc1 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 0,0,1)]+turbVisc[mapds(0, 1, 0)]+turbVisc[mapds( 0, 1, 1)])+1. / parameters.flow.Re;
  const FLOAT tmp1visc2 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 0,1,0)]+turbVisc[mapds(0, 0,-1)]+turbVisc[mapds( 0, 1,-1)])+1. / parameters.flow.Re;

  const FLOAT tmp1 = (tmp1visc1* (  (lv[mapd(0,0,1,1)] - lv[mapd(0,0,0,1)] )/lm[mapd(0,0,0,2)]
                                  + (lv[mapd(0,1,0,2)] - lv[mapd(0,0,0,2)] )/lm[mapd(0,0,0,1)] )
                     -tmp1visc2* (  (lv[mapd(0,0,0,1)] - lv[mapd(0,0,-1,1)])/lm[mapd(0,0,0,2)]
                                  + (lv[mapd(0,1,-1,2)]- lv[mapd(0,0,-1,2)])/lm[mapd(0,0,0,1)] ) )/lm[mapd(0,0,0,2)];
*/

  const FLOAT Ly1 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd(0,1,0,1)]);
  const FLOAT Lz0 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd(0,0,-1,2)]);
  const FLOAT Lz1 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd(0,0, 1,2)]);

  const FLOAT visc0 = bilinear(turbVisc[mapds( 0, 0,-1)], turbVisc[mapds( 0, 1,-1)], turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 0, 1, 0)],
                               Ly1,Lz0,0.5*lm[mapd( 0, 0,-1, 1)],0.5*lm[mapd( 0, 0,-1, 2)]) + 1.0/parameters.flow.Re;
  const FLOAT visc1 = bilinear(turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 0, 1, 0)], turbVisc[mapds( 0, 0, 1)], turbVisc[mapds( 0, 1, 1)],
                               Ly1,Lz1,0.5*lm[mapd( 0, 0, 0, 1)],0.5*lm[mapd( 0, 0, 0, 2)]) + 1.0/parameters.flow.Re;
//  if (fabs(visc0-tmp1visc2)>1.0e-12){handleError(1,"dz_dv_dz_plus_dw_dy_visc0");}
//  if (fabs(visc1-tmp1visc1)>1.0e-12){handleError(1,"dz_dv_dz_plus_dw_dy_visc1");}

  const FLOAT dv_dz0 = (lv[mapd(0,0,0,1)] - lv[mapd(0,0,-1,1)])/Lz0;
  const FLOAT dv_dz1 = (lv[mapd(0,0,1,1)] - lv[mapd(0,0, 0,1)])/Lz1;

  const FLOAT dw_dy0 = (lv[mapd(0,1,-1,2)]- lv[mapd(0,0,-1,2)])/Ly1;
  const FLOAT dw_dy1 = (lv[mapd(0,1, 0,2)]- lv[mapd(0,0, 0,2)])/Ly1;

  const FLOAT tmp2 = ( visc1*(dv_dz1+dw_dy1) - visc0*(dv_dz0+dw_dy0) )/lm[mapd(0,0,0,2)];
//  if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dz_dv_dz_plus_dw_dy");}
  return tmp2;
}


// needed for H, similar to Griebel book, Eq. (10,28); for details on implementation, see dy_du_dy_plus_dv_dx
inline FLOAT dx_dw_dx_plus_du_dz(const FLOAT * const lv, const FLOAT * const lm, const FLOAT * const turbVisc, const Parameters & parameters) {
/*  const FLOAT tmp1visc1 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 1,0,0)]+turbVisc[mapds(0, 0,1)]+turbVisc[mapds( 1, 0,1)])+1. / parameters.flow.Re;
  const FLOAT tmp1visc2 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 0,0,1)]+turbVisc[mapds(-1,0,0)]+turbVisc[mapds(-1, 0,1)])+1. / parameters.flow.Re;

  const FLOAT tmp1= (tmp1visc1* (  (lv[mapd(1,0,0,2)] - lv[mapd(0,0,0,2)] )/lm[mapd(0,0,0,0)]
                                 + (lv[mapd(0,0,1,0)] - lv[mapd(0,0,0,0)] )/lm[mapd(0,0,0,2)] )
                    -tmp1visc2* (  (lv[mapd(0,0,0,2)] - lv[mapd(-1,0,0,2)])/lm[mapd(0,0,0,0)]
                                 + (lv[mapd(-1,0,1,0)]- lv[mapd(-1,0,0,0)])/lm[mapd(0,0,0,2)] ) )/lm[mapd(0,0,0,0)];
*/

  const FLOAT Lz1 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd(0,0,1,2)]);
  const FLOAT Lx0 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd(-1,0,0,0)]);
  const FLOAT Lx1 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd( 1,0,0,0)]);

  const FLOAT visc0 = bilinear(turbVisc[mapds(-1, 0, 0)], turbVisc[mapds( 0, 0, 0)], turbVisc[mapds(-1, 0, 1)], turbVisc[mapds( 0, 0, 1)],
                               Lx0,Lz1,0.5*lm[mapd(-1, 0, 0, 0)],0.5*lm[mapd(-1, 0, 0, 2)]) + 1.0/parameters.flow.Re;
  const FLOAT visc1 = bilinear(turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 1, 0, 0)], turbVisc[mapds( 0, 0, 1)], turbVisc[mapds( 1, 0, 1)],
                               Lx1,Lz1,0.5*lm[mapd( 0, 0, 0, 0)],0.5*lm[mapd( 0, 0, 0, 2)]) + 1.0/parameters.flow.Re;
//  if (fabs(visc0-tmp1visc2)>1.0e-12){handleError(1,"dx_dw_dx_plus_du_dz_visc0");}
//  if (fabs(visc1-tmp1visc1)>1.0e-12){handleError(1,"dx_dw_dx_plus_du_dz_visc1");}

  const FLOAT dw_dx0 = (lv[mapd(0,0,0,2)] - lv[mapd(-1,0,0,2)])/Lx0;
  const FLOAT dw_dx1 = (lv[mapd(1,0,0,2)] - lv[mapd( 0,0,0,2)])/Lx1;

  const FLOAT du_dz0 = (lv[mapd(-1,0,1,0)]- lv[mapd(-1,0,0,0)])/Lz1;
  const FLOAT du_dz1 = (lv[mapd( 0,0,1,0)]- lv[mapd( 0,0,0,0)])/Lz1;

  const FLOAT tmp2 = ( visc1*(dw_dx1+du_dz1) - visc0*(dw_dx0+du_dz0) )/lm[mapd(0,0,0,0)];
//  if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dx_dw_dx_plus_du_dz");}
  return tmp2;

}


// for details on implementation, see dy_du_dy_plus_dv_dx
inline FLOAT dy_dw_dy_plus_dv_dz(const FLOAT * const lv, const FLOAT * const lm, const FLOAT * const turbVisc, const Parameters & parameters) {
/*  const FLOAT tmp1visc1 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 0,1,0)]+turbVisc[mapds(0, 0,1)]+turbVisc[mapds( 0, 1,1)])+1. / parameters.flow.Re;
  const FLOAT tmp1visc2 = 0.25*(turbVisc[mapds(0,0,0)]+turbVisc[mapds( 0,0,1)]+turbVisc[mapds(0,-1,0)]+turbVisc[mapds( 0,-1,1)])+1. / parameters.flow.Re;

  const FLOAT tmp1= (tmp1visc1* (  (lv[mapd(0,1,0,2)] - lv[mapd(0,0,0,2)] )/lm[mapd(0,0,0,1)]
                                 + (lv[mapd(0,0,1,1)] - lv[mapd(0,0,0,1)] )/lm[mapd(0,0,0,2)] )
                    -tmp1visc2* (  (lv[mapd(0,0,0,2)] - lv[mapd(0,-1,0,2)])/lm[mapd(0,0,0,1)]
                                 + (lv[mapd(0,-1,1,1)]- lv[mapd(0,-1,0,1)])/lm[mapd(0,0,0,2)] ) )/lm[mapd(0,0,0,1)];
*/

  const FLOAT Lz1 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd(0,0,1,2)]);
  const FLOAT Ly0 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd(0,-1,0,1)]);
  const FLOAT Ly1 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd(0, 1,0,1)]);

  const FLOAT visc0 = bilinear(turbVisc[mapds( 0,-1, 0)], turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 0,-1, 1)], turbVisc[mapds( 0, 0, 1)],
                               Ly0,Lz1,0.5*lm[mapd( 0,-1, 0, 1)],0.5*lm[mapd( 0,-1, 0, 2)]) + 1.0/parameters.flow.Re;
  const FLOAT visc1 = bilinear(turbVisc[mapds( 0, 0, 0)], turbVisc[mapds( 0, 1, 0)], turbVisc[mapds( 0, 0, 1)], turbVisc[mapds( 0, 1, 1)],
                               Ly1,Lz1,0.5*lm[mapd( 0, 0, 0, 1)],0.5*lm[mapd( 0, 0, 0, 2)]) + 1.0/parameters.flow.Re;
//  if (fabs(visc0-tmp1visc2)>1.0e-12){handleError(1,"dy_dw_dy_plus_dv_dz_visc0");}
//  if (fabs(visc1-tmp1visc1)>1.0e-12){handleError(1,"dy_dw_dy_plus_dv_dz_visc1");}

  const FLOAT dw_dy0 = (lv[mapd(0,0,0,2)] - lv[mapd(0,-1,0,2)])/Ly0;
  const FLOAT dw_dy1 = (lv[mapd(0,1,0,2)] - lv[mapd(0, 0,0,2)])/Ly1;

  const FLOAT dv_dz0 = (lv[mapd(0,-1,1,1)]- lv[mapd(0,-1,0,1)])/Lz1;
  const FLOAT dv_dz1 = (lv[mapd(0, 0,1,1)]- lv[mapd(0, 0,0,1)])/Lz1;

  const FLOAT tmp2 = ( visc1*(dw_dy1+dv_dz1) - visc0*(dw_dy0+dv_dz0) )/lm[mapd(0,0,0,1)];
//  if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dy_dw_dy_plus_dv_dz");}
  return tmp2; 
}



// evaluates first derivative of u w.r.t. y at the center of the grid cell.
// We therefore compute a central difference from the averaged u-values above and below
// the grid cell. This difference is divided by the distance of these values, given by
// 0.5*upper y-meshsize+0.5*lower y-meshsize + 1.0*meshsize of current cell
inline FLOAT dudy(const FLOAT * const lv, const FLOAT * const lm) {
/*	const FLOAT tmp1= 0.25
			* ((lv[mapd(-1, 1, 0, 0)] + lv[mapd(0, 1, 0, 0)])
					- (lv[mapd(-1, -1, 0, 0)] + lv[mapd(0, -1, 0, 0)]))
			/ lm[mapd(0,0,0,1)];*/
	const FLOAT tmp2= 0.5
                        * ((lv[mapd(-1, 1, 0, 0)] + lv[mapd(0, 1, 0, 0)])
                                        - (lv[mapd(-1, -1, 0, 0)] + lv[mapd(0, -1, 0, 0)]))
                        / (lm[mapd(0,0,0,1)]+0.5*lm[mapd(0,1,0,1)]+0.5*lm[mapd(0,-1,0,1)]);

//	if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dudy");}
	return tmp2;
}


// for details, see dudy
inline FLOAT dudz(const FLOAT * const lv, const FLOAT * const lm) {
/*	const FLOAT tmp1= 0.25
	       * (  (lv[mapd(-1, 0,  1, 0)] + lv[mapd( 0, 0, 1, 0)])
	          - (lv[mapd(-1, 0, -1, 0)] + lv[mapd( 0, 0,-1, 0)]) )
			/ lm[mapd(0,0,0,2)];*/
	const FLOAT tmp2= 0.5
               * (  (lv[mapd(-1, 0,  1, 0)] + lv[mapd( 0, 0, 1, 0)])
                  - (lv[mapd(-1, 0, -1, 0)] + lv[mapd( 0, 0,-1, 0)]) )
                        / (lm[mapd(0,0,0,2)]+0.5*lm[mapd(0,0,1,2)]+0.5*lm[mapd(0,0,-1,2)]);

//	if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dudz");}
	return tmp2;
}


// for details, see dudy
inline FLOAT dvdx(const FLOAT * const lv, const FLOAT * const lm) {
/*	const FLOAT tmp1= 0.25
			* ((lv[mapd( 1, 0, 0, 1)] + lv[mapd( 1,-1, 0, 1)])
			-  (lv[mapd(-1,-1, 0, 1)] + lv[mapd(-1, 0, 0, 1)]))
			/ lm[mapd(0,0,0,0)];*/
	const FLOAT tmp2= 0.5
                        * ((lv[mapd( 1, 0, 0, 1)] + lv[mapd( 1,-1, 0, 1)])
                        -  (lv[mapd(-1,-1, 0, 1)] + lv[mapd(-1, 0, 0, 1)]))
                        / (lm[mapd(0,0,0,0)] + 0.5*lm[mapd(1,0,0,0)] + 0.5*lm[mapd(-1,0,0,0)]);

//	if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dvdx");}
	return tmp2;
}


// for details, see dudy
inline FLOAT dvdz(const FLOAT * const lv, const FLOAT * const lm) {
/*  const FLOAT tmp1= 0.25
          * ((lv[mapd(0, 0, 1, 1)] + lv[mapd(0, -1, 1, 1)])
           - (lv[mapd(0,-1,-1, 1)] + lv[mapd(0,  0,-1, 1)]))
          / lm[mapd(0,0,0,2)];*/

  const FLOAT tmp2= 0.5
          * ((lv[mapd(0, 0, 1, 1)] + lv[mapd(0, -1, 1, 1)])
           - (lv[mapd(0,-1,-1, 1)] + lv[mapd(0,  0,-1, 1)]))
          / (lm[mapd(0,0,0,2)] + 0.5*lm[mapd(0,0,1,2)] + 0.5*lm[mapd(0,0,-1,2)]);

//  if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dvdz");}
  return tmp2;
}


// for details, see dudy
inline FLOAT dwdx(const FLOAT * const lv, const FLOAT * const lm) {
/*  const FLOAT tmp1= 0.25
          * ((lv[mapd( 1, 0,-1, 2)] + lv[mapd( 1, 0, 0, 2)])
           - (lv[mapd(-1, 0,-1, 2)] + lv[mapd(-1, 0, 0, 2)]))
          / lm[mapd(0,0,0,0)];*/

  const FLOAT tmp2= 0.5
          * ((lv[mapd( 1, 0,-1, 2)] + lv[mapd( 1, 0, 0, 2)])
           - (lv[mapd(-1, 0,-1, 2)] + lv[mapd(-1, 0, 0, 2)]))
          / (lm[mapd(0,0,0,0)] + 0.5*lm[mapd(1,0,0,0)] + 0.5*lm[mapd(-1,0,0,0)]);

//  if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dwdx");}
  return tmp2;
}


// for details, see dudy
inline FLOAT dwdy(const FLOAT * const lv, const FLOAT * const lm) {
/*  const FLOAT tmp1= 0.25
          * ((lv[mapd( 0, 1,-1, 2)] + lv[mapd( 0, 1, 0, 2)])
           - (lv[mapd( 0,-1,-1, 2)] + lv[mapd( 0,-1, 0, 2)]))
          / lm[mapd(0,0,0,1)];*/

  const FLOAT tmp2= 0.5
          * ((lv[mapd( 0, 1,-1, 2)] + lv[mapd( 0, 1, 0, 2)])
           - (lv[mapd( 0,-1,-1, 2)] + lv[mapd( 0,-1, 0, 2)]))
          / (lm[mapd(0,0,0,1)] + 0.5*lm[mapd(0,1,0,1)] + 0.5*lm[mapd(0,-1,0,1)]);

//  if (fabs(tmp1-tmp2)>1.0e-12){handleError(1,"dwdy");}
  return tmp2;
}


// computes the stress tensor sum in the cell midpoint (2D)
inline FLOAT StressTensorSum2D(const FLOAT * const lv, const FLOAT * const lm) {

	FLOAT S[2][2];
	S[0][0] = dudx(lv, lm);
	S[1][1] = dvdy(lv, lm);
	S[0][1] = 0.5*(dvdx(lv, lm) + dudy(lv, lm));
	S[1][0] = S[0][1];

	return S[0][0] * S[0][0] + S[1][1] * S[1][1] + 2.0 * (S[0][1] * S[0][1]);
}


// computes the stress tensor sum (3D)
inline FLOAT StressTensorSum3D(const FLOAT * const lv, const FLOAT * const lm) {
	FLOAT S[3][3];
	S[0][0] = dudx(lv, lm);
	S[1][1] = dvdy(lv, lm);
	S[2][2] = dwdz(lv, lm);

	S[0][1] = 0.5*(dvdx(lv, lm) + dudy(lv, lm));
	S[1][0] = S[0][1];

	S[1][2] = 0.5*(dvdz(lv, lm) + dwdy(lv, lm));
	S[2][1] = S[1][2];

	S[0][2] = 0.5*(dudz(lv, lm) + dwdx(lv, lm));
	S[2][0] = S[0][2];
	return S[0][0] * S[0][0] + S[1][1] * S[1][1] + S[2][2] * S[2][2]
			+ 2.0 * (S[0][1] * S[0][1] + S[0][2] * S[0][2] + S[1][2] * S[1][2]);
}


inline FLOAT computeTurbulentF2D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize,
const Parameters & parameters, const FLOAT * const turbulentViscosity, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 0)]
          + dt * (  2.0*dx_du_dx(localVelocity,localMeshsize,turbulentViscosity,parameters) + dy_du_dy_plus_dv_dx(localVelocity,localMeshsize,turbulentViscosity,parameters)
                  - du2dx(localVelocity, parameters,localMeshsize)
                  - duvdy(localVelocity, parameters,localMeshsize)
                  + parameters.environment.gx);
}

inline FLOAT computeTurbulentG2D(const FLOAT * const localVelocity,const FLOAT * const localMeshsize,
const Parameters & parameters, const FLOAT * const turbulentViscosity, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 1)]
           + dt * (  2.0*dy_dv_dy(localVelocity,localMeshsize,turbulentViscosity,parameters) + dx_dv_dx_plus_du_dy(localVelocity,localMeshsize,turbulentViscosity,parameters)
                   - duvdx(localVelocity, parameters,localMeshsize)
                   - dv2dy(localVelocity, parameters,localMeshsize)
                   + parameters.environment.gy);
}

inline FLOAT computeTurbulentF3D(const FLOAT * const localVelocity,const FLOAT * const localMeshsize,
const Parameters & parameters, const FLOAT * const turbulentViscosity, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 0)]
           + dt * (  2.0*dx_du_dx(localVelocity,localMeshsize,turbulentViscosity,parameters) + dy_du_dy_plus_dv_dx(localVelocity,localMeshsize,turbulentViscosity,parameters) + dz_du_dz_plus_dw_dx(localVelocity,localMeshsize,turbulentViscosity,parameters)
                   - du2dx(localVelocity, parameters, localMeshsize)
                   - duvdy(localVelocity, parameters, localMeshsize)
                   - duwdz(localVelocity, parameters, localMeshsize)
                   + parameters.environment.gx);
}

inline FLOAT computeTurbulentG3D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize,
const Parameters & parameters, const FLOAT * const turbulentViscosity, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 1)]
           + dt * (  2.0*dy_dv_dy(localVelocity,localMeshsize,turbulentViscosity,parameters) + dx_dv_dx_plus_du_dy(localVelocity,localMeshsize,turbulentViscosity,parameters) + dz_dv_dz_plus_dw_dy(localVelocity,localMeshsize,turbulentViscosity,parameters)
                   - dv2dy(localVelocity, parameters, localMeshsize)
                   - duvdx(localVelocity, parameters, localMeshsize)
                   - dvwdz(localVelocity, parameters, localMeshsize)
                   + parameters.environment.gy);
}

inline FLOAT computeTurbulentH3D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize,
const Parameters & parameters, const FLOAT * const turbulentViscosity, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 2)]
           + dt * (  2.0*dz_dw_dz(localVelocity,localMeshsize,turbulentViscosity,parameters) + dx_dw_dx_plus_du_dz(localVelocity,localMeshsize,turbulentViscosity,parameters) + dy_dw_dy_plus_dv_dz(localVelocity,localMeshsize,turbulentViscosity,parameters)
                   - dw2dz(localVelocity, parameters,localMeshsize)
                   - duwdx(localVelocity, parameters,localMeshsize)
                   - dvwdy(localVelocity, parameters,localMeshsize)
                   + parameters.environment.gz);
}

inline void loadLocalViscosity2D(TurbulentFlowField & flowField,
		FLOAT * const localViscosity, int i, int j) {
	for (int row = -1; row <= 1; row++) {
		for (int column = -1; column <= 1; column++) {
			localViscosity[13 + 3 * row + column] =
					flowField.getViscosityField().getScalar(i + column,
							j + row);
		}
	}
}

inline void loadLocalViscosity3D(TurbulentFlowField & flowField,
		FLOAT * const localViscosity, int i, int j, int k) {
	for (int layer = -1; layer <= 1; layer++) {
		for (int row = -1; row <= 1; row++) {
			for (int column = -1; column <= 1; column++) {
				localViscosity[13 + 9 * layer + 3 * row + column] =
						flowField.getViscosityField().getScalar(i + column, j + row,
								k + layer);
			}
		}
	}
}



#endif /* TURBULENCESTENCILFUNCTIONS_H_ */
