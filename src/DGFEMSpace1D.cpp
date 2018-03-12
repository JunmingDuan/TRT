/**
 * @file DGFEMSpace1D_GSL.cpp
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0 for 1D TRT
 * @
 * @date 2018-02-25
 */

#include <stdio.h>
#include <stdlib.h>
#include "DGFEMSpace1D.h"

DGFEMSpace1D::DGFEMSpace1D(u_int Nx, double xl, double xr)
  : Nx(Nx), xl(xl), xr(xr) {
  mesh.resize(Nx+1);
  h = (xr - xl)/Nx;
  for(u_int i = 0; i < Nx+1; ++i) {
    mesh[i] = i*h;
  }
  mu.resize(M);
  wgt.resize(M);
  I.resize(Nx); I1.resize(Nx); I2.resize(Nx); In.resize(Nx);
  T.resize(Nx); T1.resize(Nx); T2.resize(Nx); Tn.resize(Nx);
  BD_L.resize(M);
  BD_R.resize(M);
  for(u_int i = 0; i < Nx; ++i) {
    I[i].resize(M); I1[i].resize(M); I2[i].resize(M); In[i].resize(M);
    T[i].resize(K); T1[i].resize(K); T2[i].resize(K); Tn[i].resize(K);
    for(u_int m = 0; m < M; ++m) {
      I[i][m].resize(K); I1[i][m].resize(K); I2[i][m].resize(K); In[i][m].resize(K);
      BD_L[m].resize(K);
      BD_R[m].resize(K);
    }
  }
  BDRL_mat.resize(K, K);
  BDRR_mat.resize(K, K);
  BDLR_mat.resize(K, K);
  BDLL_mat.resize(K, K);
  prime_mat.resize(K, K);
  absorb_mat.resize(K, K);
  A.resize(K, K);
  rhs.resize(K);
}

void DGFEMSpace1D::BuildQuad(u_int np) {
  VEC<VEC<double> > pw;
  pw = QUADINFO[0].LGL(np);
  TemQuad.set_np(np);
  TemQuad.set_points(pw[0]);
  TemQuad.set_weight(pw[1]);
  QUADINFO.resize(Nx);
  VEC<double> gv(2), lv(2);
  lv[0] = -1, lv[1] = 1;
  for(u_int i = 0; i < Nx; ++i) {
    QUADINFO[i].set_np(np);
    gv[0] = mesh[i], gv[1] = mesh[i+1];
    QUADINFO[i].set_weight(pw[1]);
    VEC<double> x(np);
    for(u_int j = 0; j < np; ++j) {
      local_to_global(pw[0][j], lv, gv, &x[j]);
    }
    QUADINFO[i].set_points(x);
    QUADINFO[i].set_jacobi( local_to_global_jacobian(lv, gv),
        global_to_local_jacobian(lv, gv) );
  }
  u_PP_val.resize(np);
  I_PP_val.resize(np);
  T_PP_val.resize(np);
  //direction
  sum_wm = 2;
  if(M > 1) {
    pw = TemQuad.LG(M);
    for(u_int m = 0; m < M; ++m) {
      mu[m] = pw[0][m];
      wgt[m] = pw[1][m];
    }
  }
  else {
    mu[0] = 1;
    wgt[0] = 2;
  }

}

void DGFEMSpace1D::BuildMassmat() {
  VEC<double> lv(2);
  lv[0] = -1, lv[1] = 1;
  VEC<double> Pl = Poly(lv[0]);
  VEC<double> Pr = Poly(lv[1]);
  for(u_int k = 0; k < K; ++k) {
    for(u_int j = 0; j < K; ++j) {
      BDRL_mat(k,j) = Pr[j]*Pl[k];
      BDRR_mat(k,j) = Pr[j]*Pr[k];
      BDLR_mat(k,j) = Pl[j]*Pr[k];
      BDLL_mat(k,j) = Pl[j]*Pl[k];
      absorb_mat(k,j) = 0;
      prime_mat(k,j) = 0;
    }
  }

  VEC<double> x = TemQuad.points();
  VEC<double> w = TemQuad.weight();
  VEC<double> V;
  VEC<double> VP;
  for(u_int g = 0; g < x.size(); ++g) {
    V = Poly(x[g]);
    VP = PolyG(x[g]);
    for(u_int k = 0; k < K; ++k) {
      absorb_mat(k,k) += V[k]*V[k]*w[g];
      for(u_int j = 0; j < K; ++j) {
        prime_mat(k,j) += VP[k]*V[j]*w[g];
      }
    }
  }

}

void DGFEMSpace1D::Projection(u_int cell, func I0, double t, bU& u) {
  for(u_int m = 0; m < M; ++m) {
    u[m].setZero();
  }
  VEC<double> x = TemQuad.points();
  VEC<double> p = QUADINFO[cell].points();
  VEC<double> w = QUADINFO[cell].weight();
  double U;
  VEC<double> V;
  for(u_int m = 0; m < M; ++m) {
    VEC<double> basis(K, 0);
    for(u_int g = 0; g < x.size(); ++g) {
      V = Poly(x[g]);
      U = I0(mu[m], p[g], t);
      for(u_int k = 0; k < K; ++k) {
        u[m][k] += U*V[k]*w[g];
        basis[k] += V[k]*V[k]*w[g];
      }
    }
    for(u_int k = 0; k < K; ++k) {
      u[m][k] /= basis[k];//both are divided by jacobian
    }
  }
}

void DGFEMSpace1D::Projection(u_int cell, funcT T0, double t, EVEC& u) {
  u.setZero();
  VEC<double> x = TemQuad.points();
  VEC<double> p = QUADINFO[cell].points();
  VEC<double> w = QUADINFO[cell].weight();
  double U;
  VEC<double> V;
  VEC<double> basis(K, 0);
  for(u_int g = 0; g < x.size(); ++g) {
    V = Poly(x[g]);
    U = T0(p[g], t);
    for(u_int k = 0; k < K; ++k) {
      u[k] += U*V[k]*w[g];
      basis[k] += V[k]*V[k]*w[g];
    }
  }
  for(u_int k = 0; k < K; ++k) {
    u[k] /= basis[k];//both are divided by jacobian
  }
}


EVEC DGFEMSpace1D::Composition(const SOL& I, u_int cell, double x) {
  EVEC u(M);
  u.setZero();
  VEC<double> lv(2), gv(2), V;
  lv[0] = -1, lv[1] = 1;
  gv[0] = mesh[cell], gv[1] = mesh[cell+1];
  double lp;
  global_to_local(x, lv, gv, &lp);
  V = Poly(lp);
  for(u_int m = 0; m < M; ++m) {
    for(u_int k = 0; k < K; ++k) {
      u[m] += I[cell][m][k]*V[k];
    }
  }

  return u;
}

double DGFEMSpace1D::Composition(const VEC<EVEC>& T, u_int cell, double x) {
  double u(0);
  VEC<double> lv(2), gv(2), V;
  lv[0] = -1, lv[1] = 1;
  gv[0] = mesh[cell], gv[1] = mesh[cell+1];
  double lp;
  global_to_local(x, lv, gv, &lp);
  V = Poly(lp);
  for(u_int k = 0; k < K; ++k) {
    u += T[cell][k]*V[k];
  }

  return u;
}


void DGFEMSpace1D::Pk2val(const EVEC& u, VEC<double>& val) {
  VEC<double> x = TemQuad.points(), V;
  for(u_int g = 0; g < x.size(); ++g) {
    V = Poly(x[g]);
    val[g] = 0;
    for(u_int k = 0; k < K; ++k) {
      val[g] += u[k]*V[k];
    }
  }
}

void DGFEMSpace1D::init(func I0, funcT T0) {
  for(u_int i = 0; i < Nx; ++i) {
    Projection(i, I0, 0, I[i]);
    Projection(i, T0, 0, T[i]);
  }
}

double DGFEMSpace1D::cal_dt(const SOL& I) {
  //return 5.75e-5;
  return 5.75e-6;
}

int DGFEMSpace1D::forward_one_step_unsteady(const SOL& In, const SOL& I, const VEC<EVEC>& Tn, const VEC<EVEC>& T,
    func sigma_t, func q,
    double t, double dt, double* dtt, SOL& I_new, VEC<EVEC>& T_new, func BL, func BR) {
  return 0;
}

double DGFEMSpace1D::scattering_coe(const double T, const double dt, const double st) {
  return 0;
  //double tmp1 = 2.0*a*c*st*pow(T,3);
  //double tmp2 = 1 + sum_wm*dt/Cv*tmp1;
  //return tmp1/tmp2;
}

double DGFEMSpace1D::material_coe(const double T, const double Tn,
    const double dt, const double st) {
  return 0.5*a*c*st*pow(T,4);
  //double tmp1 = -1.5*a*c*st*pow(T,4);
  //double tmp2 = Tn - sum_wm*dt/Cv*tmp1;
  //return tmp1 + 2*a*c*st*pow(T,3)/(1+sum_wm*dt/Cv*2*a*c*st*pow(T,3))*tmp2;
}

double DGFEMSpace1D::temperature_numerator(const double T, const double Tn, const double scattering_I,
    const double dt, const double st) {
  return Tn + dt/Cv*scattering_I*st + 3.0*a*c*dt/Cv*st*pow(T,4);
}

double DGFEMSpace1D::temperature_denominator(const double T,
    const double dt, const double st) {
  return 1.0 + 4.0*a*c*dt/Cv*st*pow(T,3);
}

void DGFEMSpace1D::print_DG_coe(std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  VEC<double> x = TemQuad.points();
  for(u_int i = 0; i < Nx; ++i) {
    for(u_int m = 0; m < M; ++m) {
      for(u_int k = 0; k < K; ++k) {
        os << I[i][m][k] << " ";
      }
      os << "\n";
    }
    os << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

void DGFEMSpace1D::print_solution_integral(std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  VEC<double> x = TemQuad.points();
  for(u_int i = 0; i < Nx; ++i) {
    VEC<double> w = QUADINFO[i].weight();
    VEC<double> p = QUADINFO[i].points();
    for(u_int g = 0; g < x.size(); ++g) {
      //os << p[g] << " "  << w[g] << " " << Composition(I,i,p[g]).transpose() << "\n";
      os << p[g] << " "  << w[g] << " " << Composition(T,i,p[g]) << "\n";
    }
    os << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

void DGFEMSpace1D::print_solution_average(std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  for(u_int i = 0; i < Nx; ++i) {
      os << 0.5*(mesh[i]+mesh[i+1]) << " " << I[i][0] << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

