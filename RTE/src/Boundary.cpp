#include "DGFEMSpace1D.h"

void DGFEMSpace1D::Boundary(func BL, func BR, VEC<EVEC>& BD_L, VEC<EVEC>& BD_R, const double t) {
  double U1, U2;
  VEC<double> VL, VR;
  VL = Poly(-1);
  VR = Poly(1);
  for(u_int m = 0; m < M; ++m) {
    U1 = BL(mu[m], xl, t);
    U2 = BR(mu[m], xr, t);
    for(u_int k = 0; k < K; ++k) {
      BD_L[m][k] = U1*VL[k];
      BD_R[m][k] = U2*VR[k];
    }
  }
}

