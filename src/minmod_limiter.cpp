#include "minmod_limiter.h"

double minmod(const double alpha,
    const double u1, const double u2, const double u3, const double k) {
  if(k < 0 && (u2-u1) < 0 && (u3-u2) < 0) {
    return -std::min(alpha*std::min(fabs(u2-u1),fabs(u3-u2)), fabs(k));
  }
  else if(k > 0 && (u2-u1) > 0 && (u3-u2) > 0) {
    return std::min(alpha*std::min(fabs(u2-u1),fabs(u3-u2)), fabs(k));
  }
  else return 0;
}

void minmod_limiter::run(SOL& sol, const double h)
{
  u_int Nx = sol.size();
  for(u_int m = 0; m < sol[0].size(); ++m) {
    VEC<double> S(Nx);
    for(u_int i = 0; i < Nx; ++i) {
      if(i == 0) {
        S[i] = minmod(alpha,sol[i][m][1]*h,sol[i][m][0],sol[i+1][m][0],sol[i][m][1]*h);
      }
      else if(i == Nx-1) {
        S[i] = minmod(alpha,sol[i-1][m][0],sol[i][m][0],sol[i][m][1]*h,sol[i][m][1]*h);
      }
      else {
        S[i] = minmod(alpha,sol[i-1][m][0],sol[i][m][0],sol[i+1][m][0],sol[i][m][1]*h);
      }
    }
    for(u_int i = 0; i < Nx; ++i) {
      sol[i][m][1] = S[i]/h;
    }
    //
    for(u_int i = 0; i < Nx; ++i) {
      if(i == 0) {
        S[i] = minmod(alpha,sol[i][m][1]*h,sol[i][m][0],sol[i+1][m][0],sol[i][m][1]*h);
      }
      else if(i == Nx-1) {
        S[i] = minmod(alpha,sol[i-1][m][0],sol[i][m][0],sol[i][m][1]*h,sol[i][m][1]*h);
      }
      else {
        S[i] = minmod(alpha,sol[i-1][m][0],sol[i][m][0],sol[i+1][m][0],sol[i][m][1]*h);
      }
    }
    for(u_int i = 0; i < Nx; ++i) {
      sol[i][m][1] = S[i]/h;
    }

  }
}

void minmod_limiter::run(VEC<EVEC>& sol, const double h)
{
  u_int Nx = sol.size();
  VEC<double> S(Nx);
  for(u_int i = 0; i < Nx; ++i) {
    if(i == 0) {
      S[i] = minmod(alpha,sol[i][1]*h,sol[i][0],sol[i+1][0],sol[i][1]*h);
    }
    else if(i == Nx-1) {
      S[i] = minmod(alpha,sol[i-1][0],sol[i][0],sol[i][1]*h,sol[i][1]*h);
    }
    else {
      S[i] = minmod(alpha,sol[i-1][0],sol[i][0],sol[i+1][0],sol[i][1]*h);
    }
  }
  for(u_int i = 0; i < Nx; ++i) {
    sol[i][1] = S[i]/h;
  }
}

