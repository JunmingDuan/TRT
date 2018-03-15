#include "para.h"

u_int alpha = 0;
u_int M = 8;
u_int K = 5;
//0 for ghost = 0, 1 for flux = 0, 2 for period BD
u_int BDL = 1; u_int BDR = 1;
u_int PP_limiter = 0;
double EPS = 0;

double a = 0.01372;
double c = 29.98;
double Cv = 0.3;

double tol_I = 1e-14;
double tol_T = 1e-14;
double TOL = 1e-14;
u_int MAXITE = 1e4;

