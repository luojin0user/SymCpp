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

#include "EnumTypes.h"

using SymEngine::add;
using SymEngine::Basic;
using SymEngine::diff;
using SymEngine::div;
using SymEngine::function_symbol; // 用于表示 symfun/diff/int
using SymEngine::integer;
using SymEngine::max;
using SymEngine::min;
using SymEngine::mul;
using SymEngine::RCP;
using SymEngine::sin;
using SymEngine::sub;
using SymEngine::symbol;
using SymEngine::Symbol;

Case1::Case1(Region_Consts &c)
    : BasicCase(c)
{
    // obj.num_coeffs = 4;
    this->num_coeffs = 4;

    this->coeffs_exists = {0, !const_vals.Bn, 0, !const_vals.Tn, !const_vals.Ln, !const_vals.Rn};
}

// 适配父类接口
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
}

void Case1::gen_coefficient_func()
{
    // 首先需要知道传入的边界方程，然后直接带入这个方程，这样子就是直接将边界方程带入下面的等式中
    // 也就是经过这个计算之后就直接确定了最后的参数方程
    // 计算 c_hx, d_hx, e_ny, f_ny 的积分表达式
    // 这里需要注意，传入的方程 e.g. obj.B_funcs 是一个cell数组，其中的每一项是这个边界的另一边的函数的c d e f系数，
    // 对应 obj.B_coeffs ，所以需要对这个边界的4个系数分别计算，也就最多会产生16个系数方程

    // --- Bottom Boundary (c_hx) ---
    if (!this->const_vals.Bn)
        this->gen_integral(B_funcs, B_ESfuncs, 1, eq_BC[1], eq_ES[1]);

    // --- Top Boundary (d_hx) ---
    if (!this->const_vals.Tn)
        this->gen_integral(T_funcs, T_ESfuncs, 3, eq_BC[3], eq_ES[3]);

    // --- Left Boundary (e_ny) ---
    if (!this->const_vals.Ln)
        this->gen_integral(L_funcs, L_ESfuncs, 4, eq_BC[4], eq_ES[4]);

    // --- Right Boundary (f_ny) ---
    if (!this->const_vals.Rn)
        this->gen_integral(R_funcs, R_ESfuncs, 5, eq_BC[5], eq_ES[5]);
}

void Case1::gen_integral(const std::vector<Boundary_Funcs> &BC_func,
                         const std::vector<SymEngine::RCP<const SymEngine::Basic>> &BC_ESfunc,
                         int BTLR, // 指示上下左右，1开始
                         std::vector<std::array<Integral_Func, 6UL>> &eq,
                         std::vector<Integral_Func> &eq_ES,
                         std::vector<std::array<Integral_Func, 6UL>> *eq_c0d0)
{

    eq.resize(BC_func.size());

    RCP<const Basic> coeff;
    RCP<const Basic> sin_term;
    RCP<const Basic> x_or_y;

    RCP<const Basic> coeff0;

    if (BTLR == 1 || BTLR == 3)
    {
        coeff = div(integer(2), this->tau_x);
        sin_term = sin(mul(this->beta_h, sub(this->x, this->xl)));
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
        sin_term = sin(mul(this->lambda_n, sub(this->y, this->yl)));
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

            RCP<const Basic> integrand = mul(BC_func[row].funcs[col], sin_term);
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
            RCP<const Basic> es_integrand = mul(BC_ESfunc[row_ES++], sin_term);
            RCP<const Basic> es_expr = mul(coeff, es_integrand);
            // 这里可能是多个不同的积分项相加，因此分别储存所有项目，最后计算的时候再相加
            eq_ES.push_back(Integral_Func(edge_idx, es_expr, x_or_y, int_start, int_end));
        }
    }
}