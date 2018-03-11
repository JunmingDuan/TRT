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
  I.resize(Nx); I1.resize(Nx); In.resize(Nx);
  T.resize(Nx); T1.resize(Nx); Tn.resize(Nx);
  BD_L.resize(M);
  BD_R.resize(M);
  for(u_int i = 0; i < Nx; ++i) {
    I[i].resize(M); I1[i].resize(M); In[i].resize(M);
    T[i].resize(K); T1[i].resize(K); Tn[i].resize(K);
    for(u_int m = 0; m < M; ++m) {
      I[i][m].resize(K);
      I1[i][m].resize(K);
      In[i][m].resize(K);
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
  //std::cout << BDRL_mat << std::endl;
  //std::cout << BDRR_mat << std::endl;
  //std::cout << BDLR_mat << std::endl;
  //std::cout << BDLL_mat << std::endl;
  //std::cout << prime_mat << std::endl;
  //std::cout << absorb_mat << std::endl;

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
  //for(u_int i = 0; i < Nx; ++i) {
    //std::cout << I[i][0].transpose() << std::endl;
  //}
  //for(u_int i = 0; i < Nx; ++i) {
    //std::cout << T[i].transpose() << std::endl;
  //}
}

double DGFEMSpace1D::cal_dt(const SOL& I) {
  return 5.75e-6;
  //return 5.75e-6;
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

void DGFEMSpace1D::RAD_BE_unsteady(const SOL& In, const SOL& I, const VEC<EVEC>& Tn, const VEC<EVEC>& T,
    func sigma_t, func q,
    const double t, const double dt, SOL& I_new, VEC<EVEC>& T_new, func BL, func BR) {
  I_new = I;
  Boundary(BL, BR, BD_L, BD_R, t+dt);
  for(u_int m = 0; m < M; ++m) {//for each light direction
    if(mu[m] >= 0) {
      for(u_int i = 0; i < Nx; ++i) {//for each cell, update I_{m,i}^{l}
        VEC<double> x = TemQuad.points();
        VEC<double> p = QUADINFO[i].points();
        VEC<double> w = QUADINFO[i].weight();
        VEC<double> V;
        double Q;
        EVEC Im;
        A.setZero();
        rhs.setZero();
        //build matrix, outgoing and spacial prime
        A =  mu[m]*BDRR_mat
          - (mu[m])*prime_mat;//g2l_jab given by prime is eliminated by the l2g_jab given byintegral transformation
        //inflow
        if(i == 0) { rhs = mu[m]*BD_L[m]; }
        else { rhs = mu[m]*BDRL_mat*I_new[i-1][m]; }
        //time direvative
        rhs += (1./(c*dt)*QUADINFO[i].l2g_jacobian())*absorb_mat*In[i][m];
        //sigma_t*I and material, integral
        for(u_int g = 0; g < x.size(); ++g) {
          V = Poly(x[g]);
          Im = Composition(I_new, i, p[g]);
          //material temperature
          double Tg = Composition(T, i, p[g]);
          double Tng = Composition(Tn, i, p[g]);
          double st = sigma_t(mu[m],Tg,p[g]);
          for(u_int k = 0; k < K; ++k) {
            rhs[k] += material_coe(Tg, Tng, dt, st)
              *V[k]*w[g]*QUADINFO[i].l2g_jacobian();
          }
          //scattering
          for(u_int wm = 0; wm < wgt.size(); ++wm) {
            for(u_int k = 0; k < K; ++k) {
              rhs[k] += scattering_coe(Tg, dt, st)*dt/Cv*st
                *(wgt[wm]*Im[wm])
                *V[k]*w[g]*QUADINFO[i].l2g_jacobian();
            }
          }
          //I^{n+1} and sigma_t*I
          Q = 1./(c*dt)+st;
          for(u_int k = 0; k < K; ++k) {
            for(u_int j = 0; j < K; ++j) {
              A(k,j) += Q*V[j]*V[k]*w[g]*QUADINFO[i].l2g_jacobian();
            }
          }
        }

        solve_leqn(A, rhs, I_new[i][m]);
        if(PP_limiter == 1) {//perform PP limiter
        //if(PP_limiter == 2) {//perform PP limiter
          Pk2val(I_new[i][m], I_PP_val);
          u_int g = 0;
          while(g < x.size()) {
            if(I_PP_val[g] < EPS) {
              //std::cout << I_PP_val << std::endl;
              //std::cout << "i: " << i << " ,m: " << m << "\n";
              //std::cout << "before I's limiter:\n";
              //std::cout << I_new[i][m].transpose() << std::endl;
              //std::cout << "Use PP limiter on I!!" << std::endl;
              scaling_limiter::run(I_PP_val, I_new[i][m]);
              //std::cout << "after I's limiter:" << std::endl;
              //std::cout << I_new[i][m].transpose() << std::endl;
              //abort();
              break;
            }
            g++;
          }
        }

        //std::cout << "############  I  ###########" << std::endl;
        //std::cout << "i: " << i << " \nI: " << Composition(I_new,i,mesh[i]).transpose()
          //<< " \nI: " << Composition(I_new,i,mesh[i+1]).transpose() << std::endl;
        //std::cout << "I: " << Composition(I_new,i,0.5*(mesh[i]+mesh[i+1])).transpose() << std::endl;
      }
    }
    else {
      for(int i = Nx-1; i >= 0; --i) {//for each cell
        VEC<double> x = TemQuad.points();
        VEC<double> p = QUADINFO[i].points();
        VEC<double> w = QUADINFO[i].weight();
        VEC<double> V;
        double Q;
        EVEC Im;
        //build matrix, outgoing and spacial prime
        A = -mu[m]*BDLL_mat
          - (mu[m])*prime_mat;//g2l_jab given by prime is eliminated by the l2g_jab given byintegral transformation
        //inflow
        //if(i == Nx-1) {//reflex_BD
          //VEC<double> VR;
          //VR = Poly(1);
          //EVEC reflex_BD = Composition(I_new, i, mesh[i+1]);
          //for(u_int k = 0; k < K; ++k) {
            //BD_R[m][k] = reflex_BD[M-1-m]*VR[k];
          //}
          //rhs = -mu[m]*BD_R[m];
        //}
        if(i == Nx-1) {//Dirichlet_BD
          rhs = -mu[m]*BD_R[m];
        }
        else { rhs = - mu[m]*BDLR_mat*I_new[i+1][m]; }
        //time direvative
        rhs += (1./(c*dt)*QUADINFO[i].l2g_jacobian())*absorb_mat*In[i][m];
        //sigma_t*I and material, integral
        for(u_int g = 0; g < x.size(); ++g) {
          V = Poly(x[g]);
          Im = Composition(I_new, i, p[g]);
          //material temperature
          double Tg = Composition(T, i, p[g]);
          double Tng = Composition(Tn, i, p[g]);
          double st = sigma_t(mu[m],Tg,p[g]);
          for(u_int k = 0; k < K; ++k) {
            rhs[k] += material_coe(Tg, Tng, dt, st)
              *V[k]*w[g]*QUADINFO[i].l2g_jacobian();
          }
          //scattering
          for(u_int wm = 0; wm < wgt.size(); ++wm) {
            for(u_int k = 0; k < K; ++k) {
              rhs[k] += scattering_coe(Tg, dt, st)*dt/Cv*st
                *(wgt[wm]*Im[wm])
                *V[k]*w[g]*QUADINFO[i].l2g_jacobian();
            }
          }
          //sigma_t*I
          Q = 1./(c*dt)+st;
          for(u_int k = 0; k < K; ++k) {
            for(u_int j = 0; j < K; ++j) {
              A(k,j) += Q*V[j]*V[k]*w[g]*QUADINFO[i].l2g_jacobian();
            }
          }
        }

        solve_leqn(A, rhs, I_new[i][m]);

        if(PP_limiter == 1) {//perform PP limiter
        //if(PP_limiter == 2) {//perform PP limiter
          Pk2val(I_new[i][m], I_PP_val);
          u_int g = 0;
          while(g < x.size()) {
            if(I_PP_val[g] < EPS) {
              //std::cout << "i: " << i << " ,m: " << m << "\n";
              //std::cout << "before I's limiter:\n";
              //std::cout << I_new[i][m].transpose() << std::endl;
              //std::cout << "Use PP limiter on I!!" << std::endl;
              scaling_limiter::run(I_PP_val, I_new[i][m]);
              //std::cout << "after I's limiter:" << std::endl;
              //std::cout << I_new[i][m].transpose() << std::endl;
              break;
            }
            g++;
          }
        }

        //std::cout << "############  I  ###########" << std::endl;
        //std::cout << "i: " << i << " \nI: " << Composition(I_new,i,mesh[i]).transpose()
          //<< " \nI: " << Composition(I_new,i,mesh[i+1]).transpose() << std::endl;
        //std::cout << "I: " << Composition(I_new,i,0.5*(mesh[i]+mesh[i+1])).transpose() << std::endl;
      }
    }

  }
}

void DGFEMSpace1D::solve_leqn(const EMAT& A, const EVEC& rhs, EVEC& u) {
  solver.compute(A);
  u = solver.solve(rhs);
}

void DGFEMSpace1D::temperature(const SOL& I_new, const VEC<EVEC>& Tn,
    const VEC<EVEC>& T, const double dt, func sigma_t, VEC<EVEC>& T_new) {
  for(u_int i = 0; i < Nx; ++i) {
    VEC<double> x = TemQuad.points();
    VEC<double> p = QUADINFO[i].points();
    VEC<double> w = QUADINFO[i].weight();
    VEC<double> V;
    A.setZero();
    rhs.setZero();
    for(u_int g = 0; g < x.size(); ++g) {
      V = Poly(x[g]);
      double Tg = Composition(T, i, p[g]);
      double Tng = Composition(Tn, i, p[g]);
      double st = sigma_t(0,Tg,p[g]);
      EVEC Im = Composition(I_new, i, p[g]);
      double scattering_I = 0;
      for(u_int wm = 0; wm < M; ++wm) {
        scattering_I += wgt[wm] * Im[wm];
      }
      for(u_int k = 0; k < K; ++k) {
        rhs[k] += temperature_numerator(Tg, Tng, scattering_I, dt, st)
          *V[k]*w[g]*QUADINFO[i].l2g_jacobian();
        for(u_int j = 0; j < K; ++j) {
          A(k,j) += temperature_denominator(Tg, dt, st)*V[j]
            *V[k]*w[g]*QUADINFO[i].l2g_jacobian();
        }
        std::cout << "scattering_I: " << scattering_I << std::endl;
        std::cout << "temperature_numerator: " << temperature_numerator(Tg, Tng, scattering_I, dt, st) << std::endl;
        std::cout << "temperature_denominator: " << temperature_denominator(Tg, dt, st) << std::endl;
      }
    }
    solve_leqn(A, rhs, T_new[i]);
    if(PP_limiter == 1) {//perform PP limiter
      Pk2val(T_new[i], T_PP_val);
      u_int g = 0;
      //std::cout << T_PP_val << std::endl;
      while(g < x.size()) {
        if(T_PP_val[g] < EPS) {
          //std::cout << "i: " << i << "\n";
          //std::cout << "before T's limiter:\n";
          //std::cout << T_new[i].transpose() << std::endl;
          //std::cout << "Use PP limiter on T!!" << std::endl;
          scaling_limiter::run(T_PP_val, T_new[i]);
          //std::cout << "after T's limiter:" << std::endl;
          //std::cout << T_new[i].transpose() << std::endl;
          //abort();
          break;
        }
        g++;
      }
    }


    //std::cout << "############  T  ###########" << std::endl;
    //std::cout << "i: " << i << " ,T: " << Composition(T_new,i,mesh[i])
    //<< " ,T: " << Composition(T_new,i,mesh[i+1]) << std::endl;
    //std::cout << T_new[i][0] << std::endl;
    //if( Composition(T_new,i,mesh[i]) < 0 ) {
      //std::cout << "Negative temperature!" << std::endl;
      //std::cout << "i,T_new[i]: " << "\n";
      //std::cout << T_new[i] << std::endl;
      //abort();
    //}
    //std::cout << "############Tn[0],T[0],I_new[0]###########" << std::endl;
    //std::cout << Tn[0] << " " << T[0] << " " << I_new[0] << std::endl;
    //std::cout << T_new[0] << std::endl;
    //abort();
  }
}

void DGFEMSpace1D::run_unsteady(func sigma_t, func q, func BL, func BR, double t_end) {
  int ite(0), ite_T(0), ite_I(0), is_pp(0), circle(0);
  double t(0), dt(0), dtt(0), res(1), res_I(1);
  In = I; Tn = T;
  //std::cout << "Conservation: " << 2./c*Composition(I,0,0.5*(mesh[0]+mesh[1]))[0]
    //+Cv*Composition(T,0, 0.5*(mesh[0]+mesh[1]))<< std::endl;
  VEC<double> Lt, mT;
  int MAXITE(1e4);
  while ( t < t_end && ite < MAXITE ) {
    dt = cal_dt(I);
    dt = std::min(dt, t_end-t);
    res = 1; ite_T = 0;
    res_I = 1; ite_I = 0;
    while ( res_I > TOL && ite_I < MAXITE ) {//ISI iteration
      RAD_BE_unsteady(In, I, Tn, T, sigma_t, q, t, dt, I1, T1, BL, BR);
      minmod_limiter::run(I1, h);
      res_I = cal_norm_I(I, I1, 10);
      I = I1;
      std::cout << "ite_I: " << ++ite_I << ", I's res: ";
      std::cout << res_I << std::endl;
    }
    while ( res > tol && ite_T < MAXITE ) {//Newton iteration
      temperature(I, Tn, T, dt, sigma_t, T1);//update temperature
      minmod_limiter::run(T1, h);
      res = cal_norm_T(T, T1, 10);
      T = T1;
      std::cout << "ite_T: " << ++ite_T << ", T's res: " << res << std::endl;
    }
    t += dt;
    In = I; Tn = T;
    std::cout << "ite: " << ++ite << ", dt: " << dt << ", t: " << t << std::endl;
  }
}


//void DGFEMSpace1D::run_unsteady(func sigma_t, func q, func BL, func BR, double t_end) {
//  int ite(0), ite_T(0), ite_I(0), is_pp(0), circle(0);
//  double t(0), dt(0), dtt(0), res(1), res_I(1);
//  In = I; Tn = T;
//  //std::cout << "Conservation: " << 2./c*Composition(I,0,0.5*(mesh[0]+mesh[1]))[0]
//    //+Cv*Composition(T,0, 0.5*(mesh[0]+mesh[1]))<< std::endl;
//  VEC<double> Lt, mT;
//  int MAXITE(1e4);
//  while ( t < t_end && ite < MAXITE ) {
//    dt = cal_dt(I);
//    dt = std::min(dt, t_end-t);
//    res = 1; ite_T = 0;
//    while ( res > tol && ite_T < MAXITE ) {//one temporal step
//      res_I = 1; ite_I = 0;
//      //while ( (res_I > TOL || res > TOL) && ite_I < MAXITE ) {//ISI iteration
//      while ( res_I > TOL && ite_I < MAXITE ) {//ISI iteration
//        RAD_BE_unsteady(In, I, Tn, T, sigma_t, q, t, dt, I1, T1, BL, BR);
//        minmod_limiter::run(I1, h);
//        res_I = cal_norm_I(I, I1, 10);
//        I = I1;
//        std::cout << "ite_I: " << ++ite_I << ", I's res: ";
//        std::cout << res_I << std::endl;
//        //temperature(I, Tn, T, dt, sigma_t, T1);//update temperature
//        //minmod_limiter::run(T1, h);
//        //res = cal_norm_T(T, T1, 10);
//        //T = T1;
//      }
//      temperature(I, Tn, T, dt, sigma_t, T1);//update temperature
//      minmod_limiter::run(T1, h);
//      res = cal_norm_T(T, T1, 10);
//      T = T1;
//      std::cout << "ite_T: " << ++ite_T << ", T's res: " << res << std::endl;
//    }
//    t += dt;
//    In = I; Tn = T;
//    std::cout << "ite: " << ++ite << ", dt: " << dt << ", t: " << t << std::endl;
//    //Lt.push_back(t);
//    //mT.push_back(Composition(T, 0, 0.5*(mesh[0]+mesh[1])));
//    //std::cout << "Conservation: " << 2./c*Composition(I,0,0.5*(mesh[0]+mesh[1]))[0]
//      //+Cv*Composition(T,0, 0.5*(mesh[0]+mesh[1]))<< std::endl;
//  }
//  //std::cout << "err1, err2, errf:\n"
//    //<< cal_err(I, 1, t_end) << "\t"
//    //<< cal_err(I, 2, t_end) << "\t"
//    //<< cal_err(I, 0, t_end) << std::endl;
//  //for(u_int i = 0; i < Nx; ++i) {
//    //std::cout << Composition(T, i, 0.5*(mesh[i]+mesh[i+1])) << "\t";
//  //}
//  //std::cout << std::endl;
//  //Infinite medium problem
//  //std::ofstream out("t-T");
//  //out.precision(16);
//  //out << std::showpos;
//  //out.setf(std::ios::scientific);
//  //for(u_int j = 0; j < Lt.size(); ++j) {
//    //out << Lt[j] << " " << mT[j] << "\n";
//  //}
//  //out << std::endl;
//  //out.close();
//
//}


//int DGFEMSpace1D::judge_positivity(const SOL& I) {
  //for(u_int i = 0; i < Nx; ++i) {
    //cell_average is the first moment
    //cell_average[i] = I[i][0];
    //if(cell_average[i][0] < EPS) {
      //std::cout << "av[i] < 0: " << i << " " << cell_average[i][0] << std::endl;
      //return 0;
    //}
    //if(M == 3) {
      //double e = cell_average[i][2]-0.5*pow(cell_average[i][1],2)/cell_average[i][0];
      //if(e < EPS) {
        //std::cout << "e[i] < 0: " << i << " " << e << std::endl;
        //return 0;
      //}
    //}
  //}
  //return 1;
//}

double DGFEMSpace1D::cal_norm_I(const SOL& s1, const SOL& s2, int n) {
  double RES(0);
  VEC<double> norm(M,0);
  EVEC tmp1(M), tmp2(M);
  if(n == 10) {//max norm
    //for(u_int i = 0; i < Nx; ++i) {
      //VEC<double> p = QUADINFO[i].points();
      //VEC<double> w = QUADINFO[i].weight();
      //for(u_int g = 0; g < p.size(); ++g) {
        //tmp1 = Composition(s1,i,p[g]);
        //tmp2 = Composition(s2,i,p[g]);
        //for(u_int d = 0; d < M; ++d) {
          //norm[d] = std::max(fabs(tmp1[d]-tmp2[d]), norm[d]);
        //}
      //}
    //}
    for(u_int i = 0; i < Nx; ++i) {
      for(u_int d = 0; d < M; ++d) {
        norm[d] = std::max(fabs(s1[i][d][0]-s2[i][d][0]), norm[d]);
        RES = std::max(norm[d], RES);
        //norm[d] = std::max((s1[i][d]-s2[i][d]).norm(), norm[d]);
      }
    }
    return RES;
  }
  else {
    std::cout << "Wrong norm choice!" << std::endl;
    return RES;
  }
}

double DGFEMSpace1D::cal_norm_T(const VEC<EVEC>& T, const VEC<EVEC>& T1, int n) {
  double RES(0);
  VEC<double> norm(M,0);
  EVEC tmp1(M), tmp2(M);
  if(n == 2) {
    for(u_int i = 0; i < Nx; ++i) {
      VEC<double> p = QUADINFO[i].points();
      VEC<double> w = QUADINFO[i].weight();
      for(u_int g = 0; g < p.size(); ++g) {
        for(u_int d = 0; d < M; ++d) {
          norm[d] += pow(tmp1[d]-tmp2[d], 2) * w[g];
        }
      }
    }
    for(u_int d = 0; d < M; ++d) {
      norm[d] = sqrt(norm[d]/Nx);
    }
    return RES;
  }
  else if(n == 10) {//max norm
    //for(u_int i = 0; i < Nx; ++i) {
      //VEC<double> p = QUADINFO[i].points();
      //VEC<double> w = QUADINFO[i].weight();
      //for(u_int g = 0; g < p.size(); ++g) {
        //tmp1 = Composition(s1,i,p[g]);
        //tmp2 = Composition(s2,i,p[g]);
        //for(u_int d = 0; d < M; ++d) {
          //norm[d] = std::max(fabs(tmp1[d]-tmp2[d]), norm[d]);
        //}
      //}
    //}
    for(u_int i = 0; i < Nx; ++i) {
      RES = std::max(fabs(T[i][0]-T1[i][0]), RES);
      //norm[d] = std::max((s1[i][d]-s2[i][d]).norm(), norm[d]);
    }
    return RES;
  }
  else {
    std::cout << "Wrong norm choice!" << std::endl;
    return RES;
  }
}

double exact(const double mu, const double x, const double t) {
  return pow(mu,2)*pow(cos(M_PI*(x+t)), 4);
}

double DGFEMSpace1D::cal_err(const SOL& s1, int n, double t_end) {
  EVEC tmp1;
  double tmp2, NORM(0);
  double pnt;
  VEC<double> norm(M,0);
  u_int Ne(8);
  VEC<VEC<double>> pw = TemQuad.LG(Ne);
  VEC<VEC<double>> x(Nx);
  VEC<double> jab(Nx);
  for(u_int i = 0; i < Nx; ++i) {
    x[i].resize(Ne);
  }
  VEC<double> gv(2), lv(2);
  lv[0] = -1, lv[1] = 1;
  for(u_int i = 0; i < Nx; ++i) {
    gv[0] = mesh[i], gv[1] = mesh[i+1];
    for(u_int j = 0; j < Ne; ++j) {
      local_to_global(pw[0][j], lv, gv, &x[i][j]);
    }
    jab[i] = local_to_global_jacobian(lv, gv);
  }

  if(n == 2) {
    for(u_int i = 0; i < Nx; ++i) {
      for(u_int g = 0; g < Ne; ++g) {
        pnt = x[i][g];
        tmp1 = Composition(s1,i,pnt);
        for(u_int m = 0; m < M; ++m) {
          tmp2 = exact(mu[m], pnt,t_end);
          norm[m] += pow(tmp1[m]-tmp2, 2)*pw[1][g];
        }
      }
    }
    for(u_int m = 0; m < M; ++m) {
      NORM = std::max(NORM, sqrt(norm[m]/Nx));
    }
    return NORM;
  }
  else if(n == 1) {
    for(u_int i = 0; i < Nx; ++i) {
      for(u_int g = 0; g < Ne; ++g) {
        pnt = x[i][g];
        tmp1 = Composition(s1,i,pnt);
        for(u_int m = 0; m < M; ++m) {
          tmp2 = exact(mu[m], pnt,t_end);
          norm[m] += fabs(tmp1[m]-tmp2)*pw[1][g];
        }
      }
    }
    for(u_int m = 0; m < M; ++m) {
      NORM = std::max(NORM, norm[m]/Nx);
    }
    return NORM;
  }
  else if(n == 0) {
    for(u_int i = 0; i < Nx; ++i) {
      for(u_int g = 0; g < Ne; ++g) {
        pnt = x[i][g];
        tmp1 = Composition(s1,i,pnt);
        for(u_int m = 0; m < M; ++m) {
          tmp2 = exact(mu[m], pnt,t_end);
          norm[m] = std::max(norm[m], tmp1[m]-tmp2);
        }
      }
    }
    for(u_int m = 0; m < M; ++m) {
      NORM = std::max(NORM, norm[m]);
    }
    return NORM;
  }
  else {
    std::cout << "Wrong norm choice!" << std::endl;
    return NORM;
  }
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

