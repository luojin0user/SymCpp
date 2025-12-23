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
