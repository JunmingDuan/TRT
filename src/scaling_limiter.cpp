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
  t = (sol[0] - EPS) / (sol[0] - min_val);
  for(u_int k = 1; k < K; ++k) {
    sol[k] *= t;
  }
}

