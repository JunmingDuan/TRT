/**
 * @file DGFEMSpace1D.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-02-20
 */

#ifndef DGFEMSPACE1D_H
#define DGFEMSPACE1D_H

#include "Eigen/Dense"
#include "para.h"
#include "BasFun.h"
#include "Quadrature.h"
#include "interval_crd_trs.h"
#include "scaling_limiter.h"


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
    EMAT BDRL_mat;
    EMAT BDRR_mat;
    EMAT BDLR_mat;
    EMAT BDLL_mat;
    EMAT prime_mat;
    EMAT absorb_mat;
    SOL I, I1;//dimension: Nx*K*M
    SOL In;
    //std::vector<std::vector<double>> cell_average;//dimension: Nx*M
    //VEC<VEC<VEC<double>>> cell_val;//dimension: Nx*G*M
    BM bml, bmr;
    EMAT A;
    EVEC rhs;
    VEC<EVEC> BD_L, BD_R;
    Eigen::ColPivHouseholderQR<EMAT> solver;
    //, vec_u1, vec_u2;

  public:
    DGFEMSpace1D(u_int Nx, double xl, double xr);
    void BuildQuad(u_int np);
    void BuildMassmat();
    void Projection(u_int cell, func I0, double t, bU&);
    EVEC Composition(const SOL&, u_int cell, double x);
    void Pk2val(const SOL&, VEC<VEC<VEC<double>>>&);
    void init(func I0);
    double cal_dt(const SOL&, func);
    /**
     * @brief forward_one_step
     *
     * @param double dt
     *
     * @return 0, donnot change dt; 1, change dt to dtt
     */
    int forward_one_step(SOL&, func, func, func,
        SOL&, func, func);
    int forward_one_step_unsteady(const SOL&, const SOL&, func, func, func,
        double t, double dt, double* dtt, SOL&, func, func);
    //radiation part, backward Euler
    void RAD_steady(const SOL& I, func, func, func, SOL&, func, func);
    void RAD_BE_unsteady(const SOL& In, const SOL& I, func, func, func, const double, const double, SOL&, func, func);
    void Boundary(func BL, func BR, VEC<EVEC>& BD_L, VEC<EVEC>& BD_R, const double t);
    void solve_leqn(const EMAT&, const EVEC&, EVEC&);
    void run_steady(func, func, func, func, func);
    void run_unsteady(func, func, func, func, func, double t_end);
    //int judge_positivity(const SOL&);
    VEC<double> cal_norm(const SOL&, const SOL&, int);
    double cal_err(const SOL& s1, int n, double t_end);
    void print_DG_coe(std::ostream&);
    void print_solution_integral(std::ostream&);
    void print_solution_average(std::ostream&);
};

#endif //DGFEMSPACE1D_H

