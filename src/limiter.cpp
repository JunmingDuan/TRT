#include "DGFEMSpace1D.h"

void DGFEMSpace1D::do_scaling_limiter(EVEC& u) {
  Pk2val(u, u_PP_val);
  u_int g = 0;
  while(g < u_PP_val.size()) {
    if(u_PP_val[g] < EPS) {
      //std::cout << I_PP_val << std::endl;
      //std::cout << "i: " << i << " ,m: " << m << "\n";
      //std::cout << "before I's limiter:\n";
      //std::cout << I_new[i][m].transpose() << std::endl;
      //std::cout << "Use PP limiter on I!!" << std::endl;
      scaling_limiter::run(u_PP_val, u);
      //std::cout << "after I's limiter:" << std::endl;
      //std::cout << I_new[i][m].transpose() << std::endl;
      //abort();
      break;
    }
    g++;
  }
}

