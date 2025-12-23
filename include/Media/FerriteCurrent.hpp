#pragma once

#include "Case2.hpp"
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/add.h>
#include <symengine/mul.h>

class FerriteCurrent : public Case2
{
public:
    // 构造函数直接在头文件中实现
    FerriteCurrent(Region_Consts &c)
        : Case2(c)
    {
        // 局部使用命名空间，防止污染
        using SymEngine::add;
        using SymEngine::Basic;
        using SymEngine::cos;
        using SymEngine::cosh;
        using SymEngine::div;
        using SymEngine::integer;
        using SymEngine::mul;
        using SymEngine::RCP;
        using SymEngine::sin;
        using SymEngine::sinh;
        using SymEngine::sub;

        // --- 1. 沿 x 方向 A_z (A_zx_expr) ---
        // MATLAB 逻辑:
        // linear_part = (obj.yt - y) .* obj.c_0x + (y - obj.yl) .* obj.d_0x
        // harmonic_part = ( (c_hx ./ beta_h) .* sinh(beta_h .* (yt - y)) ./ sinh(beta_h .* tau_y) +
        //                   (d_hx ./ beta_h) .* sinh(beta_h .* (y - yl)) ./ sinh(beta_h .* tau_y) )
        // result = linear_part + harmonic_part .* cos(beta_h .* (x - xl))

        // 线性部分
        RCP<const Basic> term_c0 = mul(sub(this->yt, this->y), this->c_0x);
        RCP<const Basic> term_d0 = mul(sub(this->y, this->yl), this->d_0x);
        RCP<const Basic> linear_part = add(term_c0, term_d0);

        // 谐波部分公共分母: sinh(beta_h * tau_y)
        RCP<const Basic> sinh_denom_h = sinh(mul(this->beta_h, this->tau_y));

        // c_hx 项
        RCP<const Basic> term_ch_num = mul(this->c_hx, sinh(mul(this->beta_h, sub(this->yt, this->y))));
        RCP<const Basic> term_ch = div(term_ch_num, mul(this->beta_h, sinh_denom_h));

        // d_hx 项
        RCP<const Basic> term_dh_num = mul(this->d_hx, sinh(mul(this->beta_h, sub(this->y, this->yl))));
        RCP<const Basic> term_dh = div(term_dh_num, mul(this->beta_h, sinh_denom_h));

        // cos 项: cos(beta_h * (x - xl))
        RCP<const Basic> cos_part = cos(mul(this->beta_h, sub(this->x, this->xl)));

        // 组合 A_zx
        RCP<const Basic> harmonic_part = mul(add(term_ch, term_dh), cos_part);
        this->A_zx_expr = add(linear_part, harmonic_part);

        // --- 2. 沿 y 方向 A_z (A_zy_expr) ---
        // MATLAB 逻辑:
        // term_e = (-obj.e_ny ./ obj.lambda_n) .* cosh(obj.lambda_n .* (x - obj.xl)) ./ sinh(obj.lambda_n .* obj.tau_x)
        // term_f = (obj.f_ny ./ obj.lambda_n) .* cosh(obj.lambda_n .* (obj.xr - x)) ./ sinh(obj.lambda_n .* obj.tau_x)
        // result = (term_e + term_f) .* sin(obj.lambda_n .* (y - obj.yl))

        // 公共分母: sinh(lambda_n * tau_x)
        RCP<const Basic> sinh_denom_n = sinh(mul(this->lambda_n, this->tau_x));

        // e_ny 项 (注意负号)
        RCP<const Basic> term_en_coeff = div(mul(integer(-1), this->e_ny), this->lambda_n);
        RCP<const Basic> term_en_cosh = cosh(mul(this->lambda_n, sub(this->x, this->xl)));
        RCP<const Basic> term_en = div(mul(term_en_coeff, term_en_cosh), sinh_denom_n);

        // f_ny 项
        RCP<const Basic> term_fn_coeff = div(this->f_ny, this->lambda_n);
        RCP<const Basic> term_fn_cosh = cosh(mul(this->lambda_n, sub(this->xr, this->x)));
        RCP<const Basic> term_fn = div(mul(term_fn_coeff, term_fn_cosh), sinh_denom_n);

        // sin 项: sin(lambda_n * (y - yl))
        RCP<const Basic> sin_part = sin(mul(this->lambda_n, sub(this->y, this->yl)));

        // 组合 A_zy
        this->A_zy_expr = mul(add(term_en, term_fn), sin_part);
    }

    virtual ~FerriteCurrent() = default;
};