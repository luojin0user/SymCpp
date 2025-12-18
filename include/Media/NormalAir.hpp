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

class NormalAir : public Case1
{
public:
    // 构造函数直接在头文件中实现
    NormalAir(int idx,
              const SymEngine::RCP<const SymEngine::Symbol> &xl,
              const SymEngine::RCP<const SymEngine::Symbol> &xr,
              const SymEngine::RCP<const SymEngine::Symbol> &yl,
              const SymEngine::RCP<const SymEngine::Symbol> &yt,
              bool Ln, bool Rn, bool Tn, bool Bn,
              int H_max, int N_max,
              const SymEngine::RCP<const SymEngine::Symbol> &mu_r)
        : Case1(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max, mu_r)
    {
        // 局部使用命名空间，保持代码整洁且不污染全局
        using SymEngine::add;
        using SymEngine::Basic;
        using SymEngine::div;
        using SymEngine::integer;
        using SymEngine::mul;
        using SymEngine::RCP;
        using SymEngine::sin;
        using SymEngine::sinh;
        using SymEngine::sub;

        // --- 1. A_zx_expr ---
        // MATLAB 逻辑:
        // term1 = (obj.c_hx ./ obj.beta_h) .* sinh(obj.beta_h .* (obj.yt - y)) ./ sinh(obj.beta_h .* obj.tau_y)
        // term2 = (obj.d_hx ./ obj.beta_h) .* sinh(obj.beta_h .* (y - obj.yl)) ./ sinh(obj.beta_h .* obj.tau_y)
        // result = (term1 + term2) .* sin(obj.beta_h .* (x - obj.xl))

        // 公共分母: sinh(beta_h * tau_y)
        RCP<const Basic> sinh_denom_h = sinh(mul(this->beta_h, this->tau_y));

        // Term 1: c_hx 部分
        // coeff: c_hx / beta_h
        RCP<const Basic> t1_coeff = div(this->c_hx, this->beta_h);
        // sinh part: sinh(beta_h * (yt - y))
        RCP<const Basic> t1_sinh = sinh(mul(this->beta_h, sub(this->yt, this->y)));
        // Combined Term 1
        RCP<const Basic> term1 = div(mul(t1_coeff, t1_sinh), sinh_denom_h);

        // Term 2: d_hx 部分
        // coeff: d_hx / beta_h
        RCP<const Basic> t2_coeff = div(this->d_hx, this->beta_h);
        // sinh part: sinh(beta_h * (y - yl))
        RCP<const Basic> t2_sinh = sinh(mul(this->beta_h, sub(this->y, this->yl)));
        // Combined Term 2
        RCP<const Basic> term2 = div(mul(t2_coeff, t2_sinh), sinh_denom_h);

        // Sin Part: sin(beta_h * (x - xl))
        RCP<const Basic> sin_part_h = sin(mul(this->beta_h, sub(this->x, this->xl)));

        // 组合 A_zx
        this->A_zx_expr = mul(add(term1, term2), sin_part_h);

        // --- 2. A_zy_expr ---
        // MATLAB 逻辑:
        // term3 = (obj.e_ny ./ obj.lambda_n) .* sinh(obj.lambda_n .* (obj.xr - x)) ./ sinh(obj.lambda_n .* obj.tau_x)
        // term4 = (obj.f_ny ./ obj.lambda_n) .* sinh(obj.lambda_n .* (x - obj.xl)) ./ sinh(obj.lambda_n .* obj.tau_x)
        // result = (term3 + term4) .* sin(obj.lambda_n .* (y - obj.yl))

        // 公共分母: sinh(lambda_n * tau_x)
        RCP<const Basic> sinh_denom_n = sinh(mul(this->lambda_n, this->tau_x));

        // Term 3: e_ny 部分
        // coeff: e_ny / lambda_n
        RCP<const Basic> t3_coeff = div(this->e_ny, this->lambda_n);
        // sinh part: sinh(lambda_n * (xr - x))
        RCP<const Basic> t3_sinh = sinh(mul(this->lambda_n, sub(this->xr, this->x)));
        // Combined Term 3
        RCP<const Basic> term3 = div(mul(t3_coeff, t3_sinh), sinh_denom_n);

        // Term 4: f_ny 部分
        // coeff: f_ny / lambda_n
        RCP<const Basic> t4_coeff = div(this->f_ny, this->lambda_n);
        // sinh part: sinh(lambda_n * (x - xl))
        RCP<const Basic> t4_sinh = sinh(mul(this->lambda_n, sub(this->x, this->xl)));
        // Combined Term 4
        RCP<const Basic> term4 = div(mul(t4_coeff, t4_sinh), sinh_denom_n);

        // Sin Part: sin(lambda_n * (y - yl))
        RCP<const Basic> sin_part_n = sin(mul(this->lambda_n, sub(this->y, this->yl)));

        // 组合 A_zy
        this->A_zy_expr = mul(add(term3, term4), sin_part_n);
    }

    virtual ~NormalAir() = default;
};