/**
 * @file main.cpp
 * @brief Implicit DG for 1D pure advection
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-02-19
 */

#include <iostream>
#include <string>
#include <sstream>
#include "DGFEMSpace1D.h"

double I0(const double mu, const double x, const double t) {
  return 0;
}

double BL(const double mu, const double x, const double t) {
  //return 0;
  return mu*mu;
}

double BR(const double mu, const double x, const double t) {
  //return (1-exp(4*x))/2;
  return mu*mu;
}

double sigma_t(const double mu, const double x, const double t) {
  return 1;
}

double a(const double mu, const double x, const double t) {
  return 1;
}

double q(const double mu, const double x, const double t) {
  //return 1;
  return -4*pow(mu,3)*M_PI*sin(M_PI*x)*pow(cos(M_PI*x),3)
    + sigma_t(mu,x,t)*pow(mu,2)*pow(cos(M_PI*x), 4)
    - sigma_t(mu,x,t)*pow(cos(M_PI*x), 4)/3;
}

int main(int argc, char *argv[]) {
  if(argc != 4) {
    std::cout << "Usage: <Nx> <xl> <xr> " << std::endl;
    abort();
  }

  clock_t t1, t2;
  u_int Nx = atoi(argv[1]);
  double xl = atof(argv[2]);
  double xr = atof(argv[3]);
  std::cout << "Set up problem ..." << std::endl;
  DGFEMSpace1D Problem(Nx, xl, xr);
  std::cout << "Build quadrature info ..." << std::endl;
  Problem.BuildQuad(K+1);
  std::cout << "Build mass matrix ..." << std::endl;
  Problem.BuildMassmat();
  std::cout << "Initialize ..." << std::endl;
  Problem.init(I0);
  std::cout << "Start to solve ..." << std::endl;
  t1 = clock();
  Problem.run_steady(sigma_t, a, q, BL, BR);
  std::cout << "Finished ..." << std::endl;
  t2 = clock();
  std::stringstream s;
  s << "ex1_Nx" << Nx << "_K" << K << "_PP" << PP_limiter << ".dat";
  std::string filename(s.str());
  std::ofstream out(filename.c_str());
  std::cout << "Print solution to " << filename << "..." << std::endl;
  Problem.print_solution_integral(out);
  out.close();
  std::cout << "Time consumed: "
    << (t2-t1)/(double)CLOCKS_PER_SEC << std::endl;

  //exact solution
  //std::ofstream out1("ex1_exact.dat");
	//out1.precision(8);
	//out1 << std::showpos;
  //out1.setf(std::ios::scientific);
  //for(u_int i = 0; i < 1e3; ++i) {
    //double center = (xr-xl)/1e3*(i+0.5)+xl;
    //out1 << center << " " << exact(center, 0) << "\n";
  //}
  //out1 << std::endl;
  //out1 << std::defaultfloat;
  //out1.close();

  //Eigen::setNbThreads(4);
  //std::cout << Eigen::nbThreads() << std::endl;
  return 0;
}

