#include "Case2.hpp"
#include "EnumTypes.h"

#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/real_double.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/derivative.h>
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <string>
#include <iostream>

using SymEngine::add;
using SymEngine::Basic;
using SymEngine::cos;
using SymEngine::diff;
using SymEngine::div;
using SymEngine::function_symbol;
using SymEngine::integer;
using SymEngine::mul;
using SymEngine::pow;
using SymEngine::RCP;
using SymEngine::real_double;
using SymEngine::sin;
using SymEngine::sub;
using SymEngine::symbol;
using SymEngine::Symbol;

Case2::Case2(Region_Consts &c)
    : BasicCase(c)
{
    // obj.num_coeffs = 6;
    this->num_coeffs = 6;

    std::string suffix = "_" + std::to_string(idx);

    // obj.J_z = J_r;
    this->J_z = real_double(const_vals.J_r);

    // obj.c_0x = sym(...) (已在父类初始化，这里略过或重新覆盖)
    // obj.d_0x = sym(...)

    // 线性项
    // syms x y real
    // obj.A_y_P = -0.5 .* obj.mu_0 .* obj.mu_r .* obj.J_z .* y .* y;
    // -0.5 = -1/2
    RCP<const Basic> half_neg = real_double(-0.5);
    RCP<const Basic> y_sq = pow(this->y, integer(2));
    this->A_y_P = mul(mul(mul(mul(half_neg, this->mu_0), this->mu_r), this->J_z), y_sq);

    // obj.B_x_P = diff(obj.A_y_P, y);
    this->B_x_P = diff(this->A_y_P, this->y);

    // obj.B_y_P = -diff(obj.A_y_P, x);
    this->B_y_P = mul(integer(-1), diff(this->A_y_P, this->x));

    this->coeffs_exists = {!const_vals.Bn, !const_vals.Bn, !const_vals.Tn, !const_vals.Tn, !const_vals.Ln, !const_vals.Rn};
}

void Case2::gen_solution_func()
{
    // A_z = obj.A_zx_expr + obj.A_zy_expr;
    RCP<const Basic> A_z = add(this->A_zx_expr, this->A_zy_expr);

    this->B_x_x = diff(this->A_zx_expr, this->y);
    this->B_x_y = diff(this->A_zy_expr, this->y);
    this->B_y_x = mul(integer(-1), diff(this->A_zx_expr, this->x));
    this->B_y_y = mul(integer(-1), diff(this->A_zy_expr, this->x));

    RCP<const Basic> B_xx = add(this->B_x_x, this->B_x_y);
    RCP<const Basic> B_xy = add(this->B_y_x, this->B_y_y);
}

void Case2::gen_coefficient_func()
{
    // --- Bottom Boundary (c_hx) ---
    if (!this->const_vals.Bn)
        this->gen_integral(B_funcs, B_ESfuncs, 1, eq_BC[1], eq_ES[1], &eq_BC[0]);

    // --- Top Boundary (d_hx) ---
    if (!this->const_vals.Tn)
        this->gen_integral(T_funcs, T_ESfuncs, 3, eq_BC[3], eq_ES[3], &eq_BC[2]);

    // Case2的e f边界互换，换的是次序而非方程
    // --- Left Boundary (e_ny) ---
    if (!this->const_vals.Ln)
        this->gen_integral(L_funcs, L_ESfuncs, 5, eq_BC[5], eq_ES[5]);

    // --- Right Boundary (f_ny) ---
    if (!this->const_vals.Rn)
        this->gen_integral(R_funcs, R_ESfuncs, 4, eq_BC[4], eq_ES[4]);
}
