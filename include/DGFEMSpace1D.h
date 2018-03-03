/**
 * @file DGFEMSpace1D.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-02-25
 */

#ifndef DGFEMSPACE1D_H
#define DGFEMSPACE1D_H

#include "Eigen/Dense"
#include "para.h"
#include "BasFun.h"
#include "Quadrature.h"
#include "interval_crd_trs.h"
#include "scaling_limiter.h"
#include "minmod_limiter.h"


class DGFEMSpace1D {
  private:
    u_int Nx;
    double xl, xr;
    double h;
    VEC<double> mesh;
    TemplateQuadrature TemQuad;
    VEC<Quadrature> QUADINFO;
    VEC<double> mu;
    VEC<double> wgt;
    double sum_wm;
    EMAT BDRL_mat;
    EMAT BDRR_mat;
    EMAT BDLR_mat;
    EMAT BDLL_mat;
    EMAT prime_mat;
    EMAT absorb_mat;
    SOL I, I1;//dimension: Nx*K*M
    SOL In;
    VEC<EVEC> T, T1;//dimension: Nx*K
    VEC<EVEC> Tn;
    //VEC<VEC<double>> I_av;//dimension: Nx*M
    VEC<double> I_PP_val;//dimension: G
    VEC<double> T_PP_val;//dimension: G
    BM bml, bmr;
    EMAT A;
    EVEC rhs;
    VEC<EVEC> BD_L, BD_R;
    Eigen::ColPivHouseholderQR<EMAT> solver;

  public:
    DGFEMSpace1D(u_int Nx, double xl, double xr);
    void BuildQuad(u_int np);
    void BuildMassmat();
    void Projection(u_int cell, func I0, double t, bU&);
    void Projection(u_int cell, funcT T0, double t, EVEC&);
    EVEC Composition(const SOL&, u_int cell, double x);
    double Composition(const VEC<EVEC>&, u_int cell, double x);
    void Pk2val(const EVEC&, VEC<double>&);
    void init(func I0, funcT T0);
    double cal_dt(const SOL&);
    double scattering_coe(const double T, const double dt, const double st);
    double material_coe(const double T, const double Tn, const double dt, const double st);
    int forward_one_step_unsteady(const SOL&, const SOL&, const VEC<EVEC>&, const VEC<EVEC>&,
        func, func, double t, double dt, double* dtt, SOL&, VEC<EVEC>&, func, func);
    //radiation part, backward Euler
    void RAD_BE_unsteady(const SOL& In, const SOL& I, const VEC<EVEC>&, const VEC<EVEC>&,
        func, func, const double, const double, SOL&, VEC<EVEC>&, func, func);
    void temperature(const SOL& I_new, const VEC<EVEC>& Tn,
        const VEC<EVEC>& T, const double dt, func, VEC<EVEC>& T_new);
double temperature_numerator(const double T, const double Tn, const double scattering_I,
    const double dt, const double st);
double temperature_denominator(const double T,
    const double dt, const double st);
    void Boundary(func BL, func BR, VEC<EVEC>& BD_L, VEC<EVEC>& BD_R, const double t);
    void solve_leqn(const EMAT&, const EVEC&, EVEC&);
    void run_unsteady(func, func, func, func, double t_end);
    //int judge_positivity(const VEC<VEC<double>>&);
    double cal_norm_I(const SOL& s1, const SOL& s2, int n);
    double cal_norm_T(const VEC<EVEC>&, const VEC<EVEC>&, int);
    double cal_err(const SOL& s1, int n, double t_end);
    void print_DG_coe(std::ostream&);
    void print_solution_integral(std::ostream&);
    void print_solution_average(std::ostream&);
};

#endif //DGFEMSPACE1D_H

