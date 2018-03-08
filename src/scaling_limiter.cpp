#include "scaling_limiter.h"

void scaling_limiter::run(const VEC<double>& val, EVEC& sol)
{
  if(sol.size() == 1) {
    std::cout << "P1, wrong scaling_limiter choice!" << std::endl;
    return;
  }
  u_int ng = val.size();
  double t, min_val(EPS);
  for(u_int j = 0; j < ng; ++j) {
    min_val = std::min(min_val, val[j]);
  }
  if(sol[0] < 0) {
    std::cout << "average < 0! " << sol[0] << std::endl;
    abort();
  }
  //if(min_val < EPS) {
    t = (sol[0] - EPS) / (sol[0] - min_val);
    //modify moments whoes order > 1
    //std::cout << "scale: " << t << std::endl;
    for(u_int k = 1; k < K; ++k) {
      sol[k] *= t;
    }
  //}
}

//void scaling_limiter::run(const VEC<VEC<VEC<double>>>& val, SOL& sol)
//{
  //u_int n = av.size();
  //double t, min_val;

  //for(u_int i = 0; i < n; ++i) {
    //u_int m = val[i].size();
    //min_val = EPS;
    //for(u_int j = 0; j < m; ++j) {
      //if(val[i][j][0] < min_val) min_val = val[i][j][0];
    //}
    //if(min_val < EPS) {
      //std::cout << "Use PP limiter on rho!! i: " << i << std::endl;
      //t = (av[i][0] - EPS) / (av[i][0] - min_val);
      //for(u_int k = 1; k < K; ++k) {
        //sol[i][k][0] *= t;
      //}
    //}
  //}
//}

