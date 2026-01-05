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

void Case2::gen_integral(const std::vector<Boundary_Funcs> &BC_func,
                         const std::vector<SymEngine::RCP<const SymEngine::Basic>> &BC_ESfunc,
                         int BTLR, // 指示上下左右，1开始
                         std::vector<std::array<Integral_Func, 6UL>> &eq,
                         std::vector<Integral_Func> &eq_ES,
                         std::vector<std::array<Integral_Func, 6UL>> *eq_c0d0)
{

    eq.resize(BC_func.size());

    RCP<const Basic> coeff;
    RCP<const Basic> sincos_term;
    RCP<const Basic> x_or_y;

    RCP<const Basic> coeff0;

    if (BTLR == 1 || BTLR == 3)
    {
        coeff = div(integer(2), this->tau_x);
        sincos_term = cos(mul(this->beta_h, sub(this->x, this->xl)));
        x_or_y = this->x;

        if (eq_c0d0 != nullptr)
        {
            eq_c0d0->resize(BC_func.size());
            // 还需要考虑c0或者d0
            coeff0 = mul(mul(div(integer(1), this->tau_x), div(integer(1), this->tau_y)), div(integer(1), this->beta_h));
        }
    }
    else
    {
        coeff = div(integer(2), this->tau_y);
        sincos_term = sin(mul(this->lambda_n, sub(this->y, this->yl)));
        x_or_y = this->y;
    }

    for (size_t row = 0; row < BC_func.size(); row++)
    {
        size_t row_ES = 0;

        // row是对每个邻接区域的分段函数，col是对这个邻接区域的的c0 c d0 d e f
        // 对于分段函数的积分上下限，是对应的邻接区域的上下限，而不是当前区域的上下限
        // 首先找到对应临界区域
        const auto &rect = BC_func[row].rect;
        const int edge_idx = BC_func[row].idx;

        float int_start;
        float int_end;

        // 积分限
        if (BTLR == 1 || BTLR == 3)
        {
            int_start = std::max(rect.xl, this->const_vals.rect.xl);
            int_end = std::min(rect.xr, this->const_vals.rect.xr);
        }
        else
        {
            int_start = std::max(rect.yb, this->const_vals.rect.yb);
            int_end = std::min(rect.yt, this->const_vals.rect.yt);
        }

        for (size_t col = 0; col < BC_func[row].funcs.size(); col++)
        {
            if (BC_func[row].funcs[col].is_null())
                continue;

            RCP<const Basic> integrand = mul(BC_func[row].funcs[col], sincos_term);
            RCP<const Basic> expr = mul(coeff, integrand); // 积分项

            eq[row][col] = Integral_Func(edge_idx, expr, x_or_y, int_start, int_end);

            if (eq_c0d0 != nullptr)
            {
                expr = mul(BC_func[row].funcs[col], coeff0); // 积分项
                (*eq_c0d0)[row][col] = Integral_Func(edge_idx, expr, x_or_y, int_start, int_end);
            }
        }

        // BTLR B-c T-d L-e R-f
        eq_BC_loc[edge_idx] = {BTLR, row, eq_c0d0 != nullptr};

        // 处理这个区域的ES情况
        if (this->ES_regions[edge_idx])
        {
            RCP<const Basic> es_integrand = mul(BC_ESfunc[row_ES++], sincos_term);
            RCP<const Basic> es_expr = mul(coeff, es_integrand);
            // 这里可能是多个不同的积分项相加，因此分别储存所有项目，最后计算的时候再相加
            eq_ES.push_back(Integral_Func(edge_idx, es_expr, x_or_y, int_start, int_end));
        }
    }
}