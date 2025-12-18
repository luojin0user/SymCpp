#include "Case1.hpp"

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/real_double.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/derivative.h>
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <string>
#include <cmath>
#include <iostream>

using SymEngine::add;
using SymEngine::Basic;
using SymEngine::diff;
using SymEngine::div;
using SymEngine::function_symbol; // 用于表示 symfun/diff/int
using SymEngine::integer;
using SymEngine::integrate;
using SymEngine::max;
using SymEngine::min;
using SymEngine::mul;
using SymEngine::RCP;
using SymEngine::sin;
using SymEngine::sub;
using SymEngine::symbol;
using SymEngine::Symbol;

// 辅助函数：创建 SymEngine 的相等关系 (==)
RCP<const Basic> eq(const RCP<const Basic> &lhs, const RCP<const Basic> &rhs)
{
    // SymEngine 中没有直接的 "Equation" 对象用于计算，通常用 sub(lhs, rhs) == 0 或者自定义函数符号
    // 这里为了模拟 MATLAB 的 eq_A_z = (A == B)，我们使用 function_symbol "Eq"
    return function_symbol("Eq", {lhs, rhs});
}

// 辅助函数：模拟 MATLAB 的 int(expr, x, a, b)
// 真正的符号积分非常复杂，这里生成一个表示积分的符号对象
RCP<const Basic> integral_sym(const RCP<const Basic> &expr, const RCP<const Symbol> &var,
                              const RCP<const Basic> &lower, const RCP<const Basic> &upper)
{
    return function_symbol("int", {expr, var, lower, upper});
}

Case1::Case1(int idx,
             const RCP<const Symbol> &xl, const RCP<const Symbol> &xr,
             const RCP<const Symbol> &yl, const RCP<const Symbol> &yt,
             bool Ln, bool Rn, bool Tn, bool Bn,
             int H_max, int N_max,
             const RCP<const Symbol> &mu_r)
    : BasicCase(idx, xl, xr, yl, yt, H_max, N_max, mu_r)
{
    // obj.num_coeffs = 4;
    this->num_coeffs = 4;

    // obj.apply_boundaries(Ln, Rn, Tn, Bn);
    this->apply_boundaries(Ln, Rn, Tn, Bn);

    // obj.coeffs_exists = [0; ~Bn; 0; ~Tn; ~Ln; ~Rn];
    // C++ vector initialization
    this->coeffs_exists = {0, !Bn, 0, !Tn, !Ln, !Rn};

    // 初始化容器大小，防止访问越界 (根据逻辑预估)
    // 假设最大 6x6 或类似
    this->eq_c_hx.resize(10);
    for (auto &row : this->eq_c_hx)
        row.resize(10);

    this->eq_d_hx.resize(10);
    for (auto &row : this->eq_d_hx)
        row.resize(10);

    this->eq_e_ny.resize(20);
    this->eq_f_ny.resize(20);
    this->eq_ES.resize(10);
}

// 适配父类接口
void Case1::gen_solution_func(bool Ln, bool Rn, bool Tn, bool Bn)
{
    this->gen_solution_func();
}

void Case1::gen_solution_func()
{
    // syms x y (已在父类中定义)

    // A_z(x,y) = obj.A_zx_expr + obj.A_zy_expr;
    RCP<const Basic> A_z = add(this->A_zx_expr, this->A_zy_expr);

    // 计算磁场分量
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

    // 构造符号方程
    std::string suffix = "_" + std::to_string(this->idx);

    // obj.eq_A_z = symfun(...) == A_z
    RCP<const Basic> sym_Az = function_symbol("A_z" + suffix, {this->x, this->y});
    this->eq_A_z = eq(sym_Az, A_z);

    // obj.eq_B_x = symfun(...) == (B_xx + B_xy)
    RCP<const Basic> sym_Bx = function_symbol("B_x" + suffix, {this->x, this->y});
    this->eq_B_x = eq(sym_Bx, add(B_xx, B_xy));
}

void Case1::gen_coefficient_func()
{
    std::string suffix = "_" + std::to_string(this->idx);

    // 初始化累加器
    RCP<const Basic> es_c_expr = integer(0);
    RCP<const Basic> es_d_expr = integer(0);
    RCP<const Basic> es_e_expr = integer(0);
    RCP<const Basic> es_f_expr = integer(0);

    int func_num = 1;

    // obj.eq_c0x = cell(1,6); (在成员变量中处理)
    func_num++;

    int rowk = 0;
    int rowES = 1; // MATLAB index starts at 1

    // --- Bottom Boundary (c_hx) ---
    if (!this->Bn)
    {
        // MATLAB: for i = 1:length(obj.B_funcs)
        for (size_t i = 0; i < this->B_funcs.size(); ++i)
        {
            // MATLAB indices (1-based logic embedded in calculation)
            int i_matlab = i + 1;
            int row = std::ceil((double)i_matlab / 6.0);
            int col = (i_matlab + 5) % 6 + 1;

            if (this->B_funcs[i]->__eq__(*integer(0)))
            { // 简单判空/零
                // obj.eq_c_hx{row, col} = []; (Do nothing or clear)
                continue;
            }

            // 邻接区域索引
            int bottom_idx = this->bottoms[row - 1];       // C++ 0-based index
            auto bottom_i = this->all_regions[bottom_idx]; // shared_ptr

            // 积分限 intersection
            RCP<const Basic> int_start = max({bottom_i->xl, this->xl});
            RCP<const Basic> int_end = min({bottom_i->xr, this->xr});

            // c_hx_expr construction
            // (2 / tau_x) * int(func * sin(...))
            RCP<const Basic> coeff = div(integer(2), this->tau_x);
            RCP<const Basic> sin_term = sin(mul(this->beta_h, sub(this->x, this->xl)));
            RCP<const Basic> integrand = mul(this->B_funcs[i], sin_term);

            RCP<const Basic> integral_term = integral_sym(integrand, this->x, int_start, int_end);
            RCP<const Basic> c_hx_expr_val = mul(coeff, integral_term);

            // Store Equation: c_hx_suffix(x,y) == expr
            RCP<const Basic> lhs = function_symbol("c_hx" + suffix, {this->x, this->y});

            // 确保 vector 足够大
            if (this->eq_c_hx.size() <= row)
                this->eq_c_hx.resize(row + 1);
            if (this->eq_c_hx[row].size() <= col)
                this->eq_c_hx[row].resize(col + 1);
            this->eq_c_hx[row][col] = eq(lhs, c_hx_expr_val);

            if (row != rowk)
            {
                // obj.BCfuncs_loc_map(:, bottom_idx) = [2, func_num]
                this->BCfuncs_loc_map[bottom_idx] = {2, func_num};
                func_num++;
                rowk = row;

                // ES Region check
                if (bottom_idx < this->ES_regions.size() && this->ES_regions[bottom_idx])
                {
                    RCP<const Basic> es_integrand = mul(this->B_ESfuncs[rowES - 1], sin_term);
                    RCP<const Basic> es_expr = mul(coeff, integral_sym(es_integrand, this->x, int_start, int_end));
                    es_c_expr = add(es_c_expr, es_expr);
                    rowES++;
                }
            }
        }
    }
    else
    {
        func_num++;
    }

    // --- Top Boundary (d_hx) ---
    // obj.eq_d0x...
    func_num++;
    rowk = 0;
    rowES = 1;

    if (!this->Tn)
    {
        for (size_t i = 0; i < this->T_funcs.size(); ++i)
        {
            int i_matlab = i + 1;
            int row = std::ceil((double)i_matlab / 6.0);
            int col = (i_matlab + 5) % 6 + 1;

            if (this->T_funcs[i]->__eq__(*integer(0)))
                continue;

            int top_idx = this->tops[row - 1];
            auto top_i = this->all_regions[top_idx];

            RCP<const Basic> int_start = max({top_i->xl, this->xl});
            RCP<const Basic> int_end = min({top_i->xr, this->xr});

            RCP<const Basic> coeff = div(integer(2), this->tau_x);
            RCP<const Basic> sin_term = sin(mul(this->beta_h, sub(this->x, this->xl)));
            RCP<const Basic> integrand = mul(this->T_funcs[i], sin_term);

            RCP<const Basic> d_hx_expr_val = mul(coeff, integral_sym(integrand, this->x, int_start, int_end));

            RCP<const Basic> lhs = function_symbol("d_hx" + suffix, {this->x, this->y});

            if (this->eq_d_hx.size() <= row)
                this->eq_d_hx.resize(row + 1);
            if (this->eq_d_hx[row].size() <= col)
                this->eq_d_hx[row].resize(col + 1);
            this->eq_d_hx[row][col] = eq(lhs, d_hx_expr_val);

            if (row != rowk)
            {
                this->BCfuncs_loc_map[top_idx] = {4, func_num};
                func_num++;
                rowk = row;

                if (top_idx < this->ES_regions.size() && this->ES_regions[top_idx])
                {
                    RCP<const Basic> es_integrand = mul(this->T_ESfuncs[rowES - 1], sin_term);
                    RCP<const Basic> es_expr = mul(coeff, integral_sym(es_integrand, this->x, int_start, int_end));
                    es_d_expr = add(es_d_expr, es_expr);
                    rowES++;
                }
            }
        }
    }
    else
    {
        func_num++;
    }

    // --- Left Boundary (e_ny) ---
    if (!this->Ln)
    {
        for (size_t i = 0; i < this->L_funcs.size(); ++i)
        {
            if (this->L_funcs[i]->__eq__(*integer(0)))
                continue;

            // int(L_funcs * sin(lambda * (y - yl)), y, yl, yl+tau)
            RCP<const Basic> coeff = div(integer(2), this->tau_y);
            RCP<const Basic> sin_term = sin(mul(this->lambda_n, sub(this->y, this->yl)));
            RCP<const Basic> integrand = mul(this->L_funcs[i], sin_term);

            RCP<const Basic> e_ny_expr_val = mul(coeff, integral_sym(integrand, this->y, this->yl, add(this->yl, this->tau_y)));

            RCP<const Basic> lhs = function_symbol("e_ny" + suffix, {this->x, this->y});
            if (this->eq_e_ny.size() <= i)
                this->eq_e_ny.resize(i + 1);
            this->eq_e_ny[i] = eq(lhs, e_ny_expr_val);
        }

        int left_idx = this->lefts[0]; // MATLAB: lefts(1)
        if (left_idx < this->ES_regions.size() && this->ES_regions[left_idx])
        {
            RCP<const Basic> coeff = div(integer(2), this->tau_y);
            RCP<const Basic> sin_term = sin(mul(this->lambda_n, sub(this->y, this->yl)));
            RCP<const Basic> integrand = mul(this->L_ESfuncs[0], sin_term);

            RCP<const Basic> es_expr = mul(coeff, integral_sym(integrand, this->y, this->yl, add(this->yl, this->tau_y)));
            es_e_expr = add(es_e_expr, es_expr);
        }

        this->BCfuncs_loc_map[left_idx] = {5, func_num};
        func_num++;
    }
    else
    {
        func_num++;
    }

    // --- Right Boundary (f_ny) ---
    if (!this->Rn)
    {
        for (size_t i = 0; i < this->R_funcs.size(); ++i)
        {
            if (this->R_funcs[i]->__eq__(*integer(0)))
                continue;

            RCP<const Basic> coeff = div(integer(2), this->tau_y);
            RCP<const Basic> sin_term = sin(mul(this->lambda_n, sub(this->y, this->yl)));
            RCP<const Basic> integrand = mul(this->R_funcs[i], sin_term);

            RCP<const Basic> f_ny_expr_val = mul(coeff, integral_sym(integrand, this->y, this->yl, add(this->yl, this->tau_y)));

            RCP<const Basic> lhs = function_symbol("f_ny" + suffix, {this->x, this->y});
            if (this->eq_f_ny.size() <= i)
                this->eq_f_ny.resize(i + 1);
            this->eq_f_ny[i] = eq(lhs, f_ny_expr_val);
        }

        int right_idx = this->rights[0];
        if (right_idx < this->ES_regions.size() && this->ES_regions[right_idx])
        {
            RCP<const Basic> coeff = div(integer(2), this->tau_y);
            RCP<const Basic> sin_term = sin(mul(this->lambda_n, sub(this->y, this->yl)));
            RCP<const Basic> integrand = mul(this->R_ESfuncs[0], sin_term);

            RCP<const Basic> es_expr = mul(coeff, integral_sym(integrand, this->y, this->yl, add(this->yl, this->tau_y)));
            es_f_expr = add(es_f_expr, es_expr);
        }

        this->BCfuncs_loc_map[right_idx] = {6, func_num};
        func_num++;
    }
    else
    {
        func_num++;
    }

    // --- ES Equations ---
    // indices 0..5 map to MATLAB 1..6
    if (this->eq_ES.size() < 6)
        this->eq_ES.resize(6);

    // eq_ES{2} -> c_ES
    if (!es_c_expr->__eq__(*integer(0)))
    {
        RCP<const Basic> lhs = function_symbol("c_ES" + suffix, {this->x, this->y});
        this->eq_ES[1] = eq(lhs, es_c_expr);
    }

    // eq_ES{4} -> d_ES
    if (!es_d_expr->__eq__(*integer(0)))
    {
        RCP<const Basic> lhs = function_symbol("d_ES" + suffix, {this->x, this->y});
        this->eq_ES[3] = eq(lhs, es_d_expr);
    }

    // eq_ES{5} -> e_ES
    if (!es_e_expr->__eq__(*integer(0)))
    {
        RCP<const Basic> lhs = function_symbol("e_ES" + suffix, {this->x, this->y});
        this->eq_ES[4] = eq(lhs, es_e_expr);
    }

    // eq_ES{6} -> f_ES
    if (!es_f_expr->__eq__(*integer(0)))
    {
        RCP<const Basic> lhs = function_symbol("f_ES" + suffix, {this->x, this->y});
        this->eq_ES[5] = eq(lhs, es_f_expr);
    }
}