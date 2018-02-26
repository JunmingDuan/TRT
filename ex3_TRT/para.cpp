#include "para.h"

//number of the integral rule in light direction
u_int M = 1;
//number of basis function
u_int K = 1;
//0 for ghost = 0, 1 for flux = 0, 2 for period BD
//ex9
u_int BDL = 1; u_int BDR = 1;
//positivity preserving limiter
u_int PP_limiter = 0;
//parameters for Newton iteration
int MaxNt_ite = 1e1;
double Nt_tol = 1e-14;
double Nt_Ftol = 1e-14;
//tol for linear equation solver
//const double tol = 1e-13;
double tol = 1e-14;
//tol for steady solution solver
double TOL = 1e-10;
//eps for scaling_limiter
double EPS = 1e-13;
//Cv number
double Cv = 1;
//double Cv = 0.3e16;
double a = 1;

