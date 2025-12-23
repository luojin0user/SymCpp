#pragma once

#include "Case1.hpp"
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/add.h>
#include <symengine/mul.h>

class BTAir : public Case1
{
public:
    // 构造函数直接在头文件中实现
    BTAir(Region_Consts &c)
        : Case1(c)
    {
        // 在函数内部使用 using，避免污染全局命名空间
        using SymEngine::add;
        using SymEngine::Basic;
        using SymEngine::cosh;
        using SymEngine::div;
        using SymEngine::integer;
        using SymEngine::mul;
        using SymEngine::RCP;
        using SymEngine::sin;
        using SymEngine::sinh;
        using SymEngine::sub;

        // --- 构建 A_zx_expr ---
        // MATLAB 逻辑:
        // term1 = (-obj.c_hx / obj.beta_h) * sinh(obj.beta_h * (obj.yt - y)) / cosh(obj.beta_h * obj.tau_y)
        // term2 = (obj.d_hx / obj.beta_h) * sinh(obj.beta_h * (y - obj.yl)) / cosh(obj.beta_h * obj.tau_y)
        // A_zx = (term1 + term2) * sin(obj.beta_h * (x - obj.xl))

        // 公共分母: cosh(beta_h * tau_y)
        RCP<const Basic> cosh_denom = cosh(mul(this->beta_h, this->tau_y));

        // --- Term 1 ---
        // coeff: -c_hx / beta_h
        RCP<const Basic> t1_coeff = div(mul(integer(-1), this->c_hx), this->beta_h);
        // sinh_part: sinh(beta_h * (yt - y))
        RCP<const Basic> t1_sinh = sinh(mul(this->beta_h, sub(this->yt, this->y)));
        // Combined Term 1
        RCP<const Basic> term1 = div(mul(t1_coeff, t1_sinh), cosh_denom);

        // --- Term 2 ---
        // coeff: d_hx / beta_h
        RCP<const Basic> t2_coeff = div(this->d_hx, this->beta_h);
        // sinh_part: sinh(beta_h * (y - yl))
        RCP<const Basic> t2_sinh = sinh(mul(this->beta_h, sub(this->y, this->yl)));
        // Combined Term 2
        RCP<const Basic> term2 = div(mul(t2_coeff, t2_sinh), cosh_denom);

        // --- Sin Part ---
        // sin(beta_h * (x - xl))
        RCP<const Basic> sin_part = sin(mul(this->beta_h, sub(this->x, this->xl)));

        // Total Expression
        this->A_zx_expr = mul(add(term1, term2), sin_part);

        // --- 构建 A_zy_expr ---
        // obj.A_zy_expr = 0;
        this->A_zy_expr = integer(0);
    }

    virtual ~BTAir() = default;
};