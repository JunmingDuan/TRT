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

extern u_int alpha;
extern u_int M;
//number of basis function
extern u_int K;
//0 for ghost = 0, 1 for flux = 0, 2 for period BD
//ex9
extern u_int BDL; extern u_int BDR;
//positivity preserving limiter
extern u_int PP_limiter;
//eps for scaling_limiter
extern double EPS;

extern double a;
extern double c;
extern double Cv;

extern double tol_I;
extern double tol_T;
extern double TOL;
extern u_int MAXITE;

#endif //PARA_H_

