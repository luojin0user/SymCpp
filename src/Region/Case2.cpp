#include "Case2.hpp"

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

// 辅助函数 (同 Case1)
static RCP<const Basic> eq(const RCP<const Basic> &lhs, const RCP<const Basic> &rhs)
{
    return function_symbol("Eq", {lhs, rhs});
}

static RCP<const Basic> integral_sym(const RCP<const Basic> &expr, const RCP<const Symbol> &var,
                                     const RCP<const Basic> &lower, const RCP<const Basic> &upper)
{
    return function_symbol("int", {expr, var, lower, upper});
}

Case2::Case2(int idx,
             const RCP<const Symbol> &xl, const RCP<const Symbol> &xr,
             const RCP<const Symbol> &yl, const RCP<const Symbol> &yt,
             bool Ln, bool Rn, bool Tn, bool Bn,
             int H_max, int N_max,
             const RCP<const Symbol> &mu_r,
             const RCP<const Basic> &J_r)
    : BasicCase(idx, xl, xr, yl, yt, H_max, N_max, mu_r)
{
    // obj.num_coeffs = 6;
    this->num_coeffs = 6;

    // obj.apply_boundaries(Ln, Rn, Tn, Bn);
    this->apply_boundaries(Ln, Rn, Tn, Bn);

    std::string suffix = "_" + std::to_string(idx);

    // obj.J_z = J_r;
    this->J_z = J_r;

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

    // obj.eq_A_y_P = ...
    RCP<const Basic> lhs_AyP = function_symbol("A_y_P" + suffix, {this->x, this->y});
    this->eq_A_y_P = eq(lhs_AyP, this->A_y_P);

    // obj.eq_B_x_P = ...
    RCP<const Basic> lhs_BxP = function_symbol("B_x_P" + suffix, {this->x, this->y});
    this->eq_B_x_P = eq(lhs_BxP, this->B_x_P);

    // obj.eq_B_y_P = ...
    RCP<const Basic> lhs_ByP = function_symbol("B_y_P" + suffix, {this->x, this->y});
    this->eq_B_y_P = eq(lhs_ByP, this->B_y_P);

    // obj.coeffs_exists = [~Bn; ~Bn; ~Tn; ~Tn; ~Ln; ~Rn];
    // 注意：MATLAB中 ~Bn 是逻辑非。C++ vector<int> { !Bn, ... }
    this->coeffs_exists = {!Bn, !Bn, !Tn, !Tn, !Ln, !Rn};

    // 初始化容器
    this->eq_c_hx.resize(20);
    this->eq_d_hx.resize(20);
    this->eq_c0x.resize(20);
    this->eq_d0x.resize(20);
    this->eq_e_ny.resize(20);
    this->eq_f_ny.resize(20);
    this->eq_ES.resize(10);
}

void Case2::gen_solution_func(bool Ln, bool Rn, bool Tn, bool Bn)
{
    // A_z = obj.A_zx_expr + obj.A_zy_expr;
    RCP<const Basic> A_z = add(this->A_zx_expr, this->A_zy_expr);

    // obj.B_x_x = diff(obj.A_zx_expr , y);
    this->B_x_x = diff(this->A_zx_expr, this->y);
    // obj.B_x_y = diff(obj.A_zy_expr , y);
    this->B_x_y = diff(this->A_zy_expr, this->y);
    // obj.B_y_x = -diff(obj.A_zx_expr , x);
    this->B_y_x = mul(integer(-1), diff(this->A_zx_expr, this->x));
    // obj.B_y_y = -diff(obj.A_zy_expr , x);
    this->B_y_y = mul(integer(-1), diff(this->A_zy_expr, this->x));

    // B_xx = obj.B_x_x + obj.B_x_y;
    RCP<const Basic> B_xx = add(this->B_x_x, this->B_x_y);
    // B_xy = obj.B_y_x + obj.B_y_y;
    RCP<const Basic> B_xy = add(this->B_y_x, this->B_y_y);

    std::string suffix = "_" + std::to_string(this->idx);

    // obj.eq_A_z = symfun(...) == A_z;
    RCP<const Basic> lhs_Az = function_symbol("A_z" + suffix, {this->x, this->y});
    this->eq_A_z = eq(lhs_Az, A_z);

    // obj.eq_B_x = symfun(...) == (B_xx + B_xy);
    RCP<const Basic> lhs_Bx = function_symbol("B_x" + suffix, {this->x, this->y});
    this->eq_B_x = eq(lhs_Bx, add(B_xx, B_xy));
}

void Case2::gen_coefficient_func()
{
    std::string suffix = "_" + std::to_string(this->idx);
    int func_num = 1;

    // obj.eq_c0x = cell(1,6); ... initialized in constructor

    // --- Bottom Boundary (c_hx, c_0x) ---
    // if ~obj.Bn
    if (!this->Bn)
    {
        // es_c_expr accumulators are calculated differently in Case2 (see end of function),
        // but here we check for ES_regions to update map?
        RCP<const Basic> es_c_expr_val = integer(0); // local accumulator logic similar to Case1 if needed?
        // MATLAB Code: es_c_expr = es_c_expr + es_expr; (variable scope is function wide)

        for (size_t i = 0; i < this->B_funcs.size(); ++i)
        {
            if (this->B_funcs[i]->__eq__(*integer(0)))
            {
                // obj.eq_c_hx{i} = [];
                continue;
            }

            // c_hx_expr = (2/tau_x) * int(func * cos(beta_h * (x - xl)))
            RCP<const Basic> coeff_h = div(integer(2), this->tau_x);
            RCP<const Basic> cos_term = cos(mul(this->beta_h, sub(this->x, this->xl)));
            RCP<const Basic> integrand_h = mul(this->B_funcs[i], cos_term);

            RCP<const Basic> c_hx_int = integral_sym(integrand_h, this->x, this->xl, add(this->xl, this->tau_x));
            RCP<const Basic> c_hx_res = mul(coeff_h, c_hx_int);

            RCP<const Basic> lhs_ch = function_symbol("c_hx" + suffix, {this->x, this->y});
            if (this->eq_c_hx.size() <= i)
                this->eq_c_hx.resize(i + 1);
            this->eq_c_hx[i] = eq(lhs_ch, c_hx_res);

            // c_0x_expr = (1/tau_x) * int((1/tau_y) * (1/beta_h) * func)
            // Note: 1/beta_h is symbolic here.
            RCP<const Basic> term1 = div(integer(1), this->tau_x);
            RCP<const Basic> term2 = div(integer(1), this->tau_y);
            RCP<const Basic> term3 = div(integer(1), this->beta_h);
            RCP<const Basic> integrand_0 = mul(mul(term2, term3), this->B_funcs[i]);

            RCP<const Basic> c_0x_int = integral_sym(integrand_0, this->x, this->xl, add(this->xl, this->tau_x));
            RCP<const Basic> c_0x_res = mul(term1, c_0x_int);

            RCP<const Basic> lhs_c0 = function_symbol("c_0x" + suffix, {this->x, this->y});
            if (this->eq_c0x.size() <= i)
                this->eq_c0x.resize(i + 1);
            this->eq_c0x[i] = eq(lhs_c0, c_0x_res);
        }

        int bottom_idx = this->bottoms[0]; // MATLAB: bottoms(1)
        if (bottom_idx < this->ES_regions.size() && this->ES_regions[bottom_idx])
        {
            // es_expr = (1/tau_x) * int((1/tau_y) * B_ESfuncs(1))
            RCP<const Basic> term1 = div(integer(1), this->tau_x);
            RCP<const Basic> term2 = div(integer(1), this->tau_y);
            RCP<const Basic> integrand = mul(term2, this->B_ESfuncs[0]);

            // es_c_expr logic: MATLAB accumulates this into a variable 'es_c_expr'.
            // In C++, we need to store this.
            // Note: The MATLAB code later calculates es_c_expr again using A_y_P at the end of function.
            // Wait, look at MATLAB: `es_c_expr = 0` is NOT initialized at top of Case2.
            // But `es_c_expr = es_c_expr + es_expr` implies it exists.
            // HOWEVER, at the end of function: `es_c_expr(x,y) = ... obj.A_y_P ...` overwrites it?
            // Actually, the variable name at the end is the same. This is potentially a bug in the provided MATLAB
            // or `es_c_expr` is used for boundary accumulation, and THEN `eq_ES{1}` is set using `A_y_P`.
            // LET'S LOOK STRICTLY:
            // End of function: `es_c_expr(x,y) = (1/tau_x) * int(...) * A_y_P ...`
            // Then `obj.eq_ES{1} = ... == es_c_expr`.
            // The accumulation inside the `if ~obj.Bn` loop seems to calculate a value that is strictly
            // used for `c_ES` (index 2 in Case1) but here it's named `es_c_expr`.
            // BUT Case2 specifically sets `obj.eq_ES{2} = []`.
            // So the accumulation inside the loop might be unused code in Case2?
            // Re-reading: `es_c_expr = es_c_expr + es_expr` accumulates.
            // BUT at the end: `es_c_expr(x,y) = ...` (Assigns NEW symbolic function).
            // It effectively overwrites. I will follow the code: Implement loop logic, then implement end logic.
        }

        func_num++;
        this->BCfuncs_loc_map[bottom_idx] = {2, func_num};
        func_num++;
    }
    else
    {
        func_num++;
    }

    // --- Top Boundary (d_hx, d_0x) ---
    if (!this->Tn)
    {
        for (size_t i = 0; i < this->T_funcs.size(); ++i)
        {
            if (this->T_funcs[i]->__eq__(*integer(0)))
                continue;

            // d_hx = (2/tau_x) * int(func * cos(...))
            RCP<const Basic> coeff = div(integer(2), this->tau_x);
            RCP<const Basic> cos_term = cos(mul(this->beta_h, sub(this->x, this->xl)));
            RCP<const Basic> integrand = mul(this->T_funcs[i], cos_term);

            RCP<const Basic> d_hx_res = mul(coeff, integral_sym(integrand, this->x, this->xl, add(this->xl, this->tau_x)));

            RCP<const Basic> lhs_dh = function_symbol("d_hx" + suffix, {this->x, this->y});
            if (this->eq_d_hx.size() <= i)
                this->eq_d_hx.resize(i + 1);
            this->eq_d_hx[i] = eq(lhs_dh, d_hx_res);

            // d_0x = (1/tau_x) * int((1/tau_y)*(1/beta_h)*func)
            RCP<const Basic> term1 = div(integer(1), this->tau_x);
            RCP<const Basic> term2 = div(integer(1), this->tau_y);
            RCP<const Basic> term3 = div(integer(1), this->beta_h);
            RCP<const Basic> integrand_0 = mul(mul(term2, term3), this->T_funcs[i]);

            RCP<const Basic> d_0x_res = mul(term1, integral_sym(integrand_0, this->x, this->xl, add(this->xl, this->tau_x)));

            RCP<const Basic> lhs_d0 = function_symbol("d_0x" + suffix, {this->x, this->y});
            if (this->eq_d0x.size() <= i)
                this->eq_d0x.resize(i + 1);
            this->eq_d0x[i] = eq(lhs_d0, d_0x_res);
        }

        int top_idx = this->tops[0];
        if (top_idx < this->ES_regions.size() && this->ES_regions[top_idx])
        {
            // Accumulation logic for es_d_expr (potentially overwritten later)
        }

        func_num++;
        this->BCfuncs_loc_map[top_idx] = {4, func_num};
        func_num++;
    }
    else
    {
        func_num++;
    }

    // --- Right Boundary maps to e_ny (Left Logic in Case1, Swapped here) ---
    // if ~obj.Rn
    if (!this->Rn)
    {
        for (size_t i = 0; i < this->R_funcs.size(); ++i)
        {
            if (this->R_funcs[i]->__eq__(*integer(0)))
                continue;

            // e_ny = (2/tau_y) * int(R_funcs * sin(...))
            RCP<const Basic> coeff = div(integer(2), this->tau_y);
            RCP<const Basic> sin_term = sin(mul(this->lambda_n, sub(this->y, this->yl)));
            RCP<const Basic> integrand = mul(this->R_funcs[i], sin_term);

            RCP<const Basic> e_ny_res = mul(coeff, integral_sym(integrand, this->y, this->yl, add(this->yl, this->tau_y)));

            RCP<const Basic> lhs = function_symbol("e_ny" + suffix, {this->x, this->y});
            if (this->eq_e_ny.size() <= i)
                this->eq_e_ny.resize(i + 1);
            this->eq_e_ny[i] = eq(lhs, e_ny_res);
        }

        int right_idx = this->rights[0];
        this->BCfuncs_loc_map[right_idx] = {5, func_num};
        func_num++;
    }
    else
    {
        func_num++;
    }

    // --- Left Boundary maps to f_ny ---
    // if ~obj.Ln
    if (!this->Ln)
    {
        for (size_t i = 0; i < this->L_funcs.size(); ++i)
        {
            if (this->L_funcs[i]->__eq__(*integer(0)))
                continue;

            // f_ny = (2/tau_y) * int(L_funcs * sin(...))
            RCP<const Basic> coeff = div(integer(2), this->tau_y);
            RCP<const Basic> sin_term = sin(mul(this->lambda_n, sub(this->y, this->yl)));
            RCP<const Basic> integrand = mul(this->L_funcs[i], sin_term);

            RCP<const Basic> f_ny_res = mul(coeff, integral_sym(integrand, this->y, this->yl, add(this->yl, this->tau_y)));

            RCP<const Basic> lhs = function_symbol("f_ny" + suffix, {this->x, this->y});
            if (this->eq_f_ny.size() <= i)
                this->eq_f_ny.resize(i + 1);
            this->eq_f_ny[i] = eq(lhs, f_ny_res);
        }

        int left_idx = this->lefts[0];
        this->BCfuncs_loc_map[left_idx] = {6, func_num};
        func_num++;
    }
    else
    {
        func_num++;
    }

    // --- ES Terms (Defined at end of MATLAB function) ---
    // es_c_expr = (1/tau_x) * int((1/tau_y) * A_y_P, x, ...)
    // subs(es_c_expr, y, yl)

    RCP<const Basic> term1 = div(integer(1), this->tau_x);
    RCP<const Basic> term2 = div(integer(1), this->tau_y);
    RCP<const Basic> integrand_es = mul(term2, this->A_y_P);

    // Integral
    RCP<const Basic> es_c_int = integral_sym(integrand_es, this->x, this->xl, add(this->xl, this->tau_x));
    RCP<const Basic> es_c_val = mul(term1, es_c_int);

    // Subs y -> yl
    es_c_val = es_c_val->subs({{this->y, this->yl}});

    // obj.eq_ES{1} = ...
    RCP<const Basic> lhs_c_es = function_symbol("c_ES" + suffix, {this->x, this->y});
    this->eq_ES[0] = eq(lhs_c_es, es_c_val); // MATLAB {1} -> C++ [0]

    // es_d_expr calculation (same formula, subs y -> yt)
    RCP<const Basic> es_d_val = mul(term1, es_c_int); // Same integral expression before subs
    es_d_val = es_d_val->subs({{this->y, this->yt}});

    // obj.eq_ES{3} = ...
    RCP<const Basic> lhs_d_es = function_symbol("d_ES" + suffix, {this->x, this->y});
    this->eq_ES[2] = eq(lhs_d_es, es_d_val); // MATLAB {3} -> C++ [2]
}