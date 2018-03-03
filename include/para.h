#ifndef PARA_H_
#define PARA_H_

#include <cstdlib>
#include <vector>
#include "Eigen/Dense"
#include "VEC.h"

typedef int BM;
typedef VEC<VEC<double>> QUAD;
typedef double (*func)(const double, const double, const double);
typedef double (*funcT)(const double, const double);
typedef Eigen::MatrixXd EMAT;
typedef Eigen::VectorXd EVEC;
typedef VEC<EVEC> bU;
typedef VEC<bU> SOL;

//speed of light
extern double c;

extern u_int alpha;
extern u_int M;
//number of basis function
extern u_int K;
//0 for ghost = 0, 1 for flux = 0, 2 for period BD
//ex9
extern u_int BDL; extern u_int BDR;
//positivity preserving limiter
extern u_int PP_limiter;
//parameters for Newton iteration
extern int MaxNt_ite;
extern double Nt_tol;
extern double Nt_Ftol;
//tol for linear equation solver
//const double tol = 1e-13;
extern double tol;
//tol for steady solution solver
extern double TOL;
//eps for scaling_limiter
extern double EPS;
extern double a;
extern double Cv;

#endif //PARA_H_

