#include "DGFEMSpace1D.h"

double DGFEMSpace1D::RAD_BE_unsteady(const SOL& In, const SOL& I, const VEC<EVEC>& Tn, const VEC<EVEC>& T,
    func sigma_t, func q,
    const double t, const double dt, SOL& I_new, func BL, func BR) {
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
        //if(m == 1 && i == 0) {
          //std::cout << "A:" << std::endl;
          //std::cout << A << std::endl;
          //std::cout << "rhs:" << std::endl;
          //std::cout << rhs << std::endl;
        //}

        solve_leqn(A, rhs, I_new[i][m]);
        //std::cout << "RAD_BE_unsteady1" << std::endl;
        //std::cout << "m: " << m << std::endl;
        //std::cout << I_new << std::endl;
        if(PP_limiter == 1) {//perform PP limiter
          do_scaling_limiter(I_new[i][m]);
        }
      }//end i
    }//end ->
    else {
      u_int i = Nx-1;
      while(i >= 0) {
      //for(int i = Nx-1; i >= 0; --i) {//for each cell
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
        //std::cout << "RAD_BE_unsteady1" << std::endl;
        //std::cout << "m: " << m << std::endl;
        //std::cout << I_new << std::endl;
        if(PP_limiter == 1) {//perform PP limiter
          do_scaling_limiter(I_new[i][m]);
        }
        if(i-- == 0) break;
      }//end i
    }//end <-
  }//end M
  return cal_norm_I(I, I_new, 10);
}

void DGFEMSpace1D::solve_leqn(const EMAT& A, const EVEC& rhs, EVEC& u) {
  solver.compute(A);
  u = solver.solve(rhs);
}

double DGFEMSpace1D::temperature(const SOL& I_new, const VEC<EVEC>& Tn,
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
      }
    }
    solve_leqn(A, rhs, T_new[i]);
    //if(i == 0) {
      //std::cout << "T i: " << i << std::endl;
      //std::cout << "T A: " << std::endl;
      //std::cout << A << std::endl;
      //std::cout << "T rhs: " << std::endl;
      //std::cout << rhs << std::endl;
      //std::cout << "T_new: " << std::endl;
      //std::cout << T_new[i] << std::endl;
    //}

    if(PP_limiter == 1) {//perform PP limiter
      do_scaling_limiter(T_new[i]);
    }
  }
  return cal_norm_T(T, T_new, 10);
}

void DGFEMSpace1D::RAD_BE_unsteady_ite(const SOL& In, const SOL& I, const VEC<EVEC>& Tn, const VEC<EVEC>& T,
    func sigma_t, func q,
    const double t, const double dt, SOL& I_new, func BL, func BR) {
  double res_I = 1;
  u_int ite_I = 0;
  I2 = I;
  while ( res_I > tol_I && ite_I < MAXITE ) {//ISI iteration
    res_I = RAD_BE_unsteady(In, I2, Tn, T, sigma_t, q, t, dt, I_new, BL, BR);
    if(alpha >= 0) {
      minmod_limiter::run(I_new, h);
    }
    I2 = I_new;
    std::cout << "ite_I: " << ++ite_I << ", res_I: " << res_I << std::endl;
  }
}

void DGFEMSpace1D::temperature_ite(const SOL& I_new, const VEC<EVEC>& Tn,
    const VEC<EVEC>& T, const double dt, func sigma_t, VEC<EVEC>& T_new) {
  double res_T = 1;
  u_int ite_T = 0;
  T2 = T;
  while ( res_T > tol_T && ite_T < 1) {//Newton iteration
    res_T = temperature(I_new, Tn, T2, dt, sigma_t, T_new);
    if(alpha >= 0) {
      minmod_limiter::run(T_new, h);
    }
    T2 = T_new;
    std::cout << "ite_T: " << ++ite_T << ", res_T: " << res_T << std::endl;
  }
}

void DGFEMSpace1D::run_unsteady(func sigma_t, func q, func BL, func BR, double t_end) {
  double t(0), dt(0);
  In = I; Tn = T;
  u_int IT(0);
  //std::cout << "Conservation: " << 2./c*Composition(I,0,0.5*(mesh[0]+mesh[1]))[0]
  //+Cv*Composition(T,0, 0.5*(mesh[0]+mesh[1]))<< std::endl;
  //while (t < t_end && IT < 1) {
  while (t < t_end) {
    double res = 1;
    u_int ite= 0;
    //std::cout << "begin" << std::endl;
    //std::cout << "I" << std::endl;
    //std::cout << I << std::endl;
    //std::cout << "begin" << std::endl;
    //std::cout << "T" << std::endl;
    //std::cout << T << std::endl;
    while ( res > TOL && ite < MAXITE ) {//for a time step
      dt = cal_dt(I);
      dt = std::min(dt, t_end-t);
      RAD_BE_unsteady_ite(In, I, Tn, T, sigma_t, q, t, dt, I1, BL, BR);
      //std::cout << "ite_I" << std::endl;
      //std::cout << "I" << std::endl;
      //std::cout << I1 << std::endl;
      temperature_ite(I1, Tn, T, dt, sigma_t, T1);
      //std::cout << "one ite_T" << std::endl;
      //std::cout << "T" << std::endl;
      //std::cout << T1 << std::endl;
      res = cal_norm(I, I1, T, T1, 10);
      I = I1;
      T = T1;
      std::cout << "ite: " << ++ite << ", res: " << res << std::endl;
    }
    t += dt;
    In = I; Tn = T;
    //std::cout << In << std::endl;
    //std::cout << Tn << std::endl;
    std::cout << "t: " << t << ", dt: " << dt << std::endl;
    IT++;
  }
}

