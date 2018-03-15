/**
 * @file main.cpp
 * @brief Implicit DG for 1D Marshak wave
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-03-15
 */

#include <iostream>
#include <string>
#include <sstream>
#include "DGFEMSpace1D.h"

double T0(const double x, const double t) {
  //if(x <= 0.55/15) return 1;
  //else return 1e-1;
  //return 1e-2;
  return 1e-2;
}

double I0(const double mu, const double x, const double t) {
  return 0.5*a*c*pow(T0(x,t),4);
}

double BL(const double mu, const double x, const double t) {
  return 0.5*a*c*pow(1,4);
}

double BR(const double mu, const double x, const double t) {
  return 0.5*a*c*pow(1e-2,4);
}

double sigma_t(const double mu, const double T, const double x) {
  return 300.0/pow(T,3);
  //return 300;
  //return 10;
}

double q(const double mu, const double x, const double t) {
  return 0;
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
  Problem.init(I0, T0);
  std::cout << "Start to solve ..." << std::endl;
  t1 = clock();
  //Problem.run_unsteady(sigma_t, q, BL, BR, 0.00092);
  //Problem.run_unsteady(sigma_t, q, BL, BR, 0.0001725);
  //Problem.run_unsteady(sigma_t, q, BL, BR, 4.6e-2);
  Problem.run_unsteady(sigma_t, q, BL, BR, 2.3e-2);
  //Problem.run_unsteady(sigma_t, q, BL, BR, 5.75e-2);
  //Problem.run_unsteady(sigma_t, q, BL, BR, 1.72e-1);
  std::cout << "Finished ..." << std::endl;
  t2 = clock();
  std::stringstream s;
  s << "ex2_Nx" << Nx << "_K" << K << "_PP" << PP_limiter <<
    "_MD" << alpha << ".dat";
  std::string filename(s.str());
  std::ofstream out(filename.c_str());
  std::cout << "Print solution to " << filename << "..." << std::endl;
  //Problem.print_DG_coe(out);
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

