#include "BasicCase.hpp"

#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/real_double.h>
#include <symengine/constants.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <stdexcept>

// 在 cpp 文件中使用 using，不污染头文件
using SymEngine::Basic;
using SymEngine::div;
using SymEngine::integer;
using SymEngine::max;
using SymEngine::min;
using SymEngine::mul;
using SymEngine::pi;
using SymEngine::RCP;
using SymEngine::real_double;
using SymEngine::sub;
using SymEngine::symbol;
using SymEngine::Symbol;

SymEngine::RCP<const SymEngine::Symbol> gx = SymEngine::symbol("x");
SymEngine::RCP<const SymEngine::Symbol> gy = SymEngine::symbol("y");

BasicCase::BasicCase(Region_Consts &c) : const_vals(c)
{
    this->idx = c.idx;

    // syms x y
    this->x = gx;
    this->y = gy;

    std::string suffix = "_" + std::to_string(idx);

    // --- 公共编号变量 ---
    this->xr = real_double(const_vals.rect.xr);
    this->xl = real_double(const_vals.rect.xl);
    this->yt = real_double(const_vals.rect.yt);
    this->yl = real_double(const_vals.rect.yb);

    // obj.tau_x = obj.xr - obj.xl;
    this->tau_x = real_double(const_vals.rect.xr - const_vals.rect.xl);
    this->tau_y = real_double(const_vals.rect.yt - const_vals.rect.yb);

    // --- 公共求和指标 ---
    this->h = symbol("h" + suffix);
    this->n = symbol("n" + suffix);

    // --- 公共参数 ---
    // obj.beta_h = obj.h * pi / obj.tau_x;
    this->beta_h = div(mul(this->h, pi), this->tau_x);

    // obj.lambda_n = obj.n * pi / obj.tau_y;
    this->lambda_n = div(mul(this->n, pi), this->tau_y);

    // obj.mu_0 = 4*pi*1e-7;
    this->mu_0 = mul(mul(integer(4), pi), real_double(1e-7));
    this->mu_r = real_double(const_vals.mu_r);

    // --- 边界积分项 ---
    this->c_0x = symbol("c_0x" + suffix);
    this->d_0x = symbol("d_0x" + suffix);
    this->c_hx = symbol("c_hx" + suffix);
    this->d_hx = symbol("d_hx" + suffix);
    this->e_ny = symbol("e_ny" + suffix);
    this->f_ny = symbol("f_ny" + suffix);

    this->num_coeffs = 6;

    apply_boundaries();

    eq_BC_loc.resize(c.all_regions_num);
}

void BasicCase::apply_boundaries()
{
    if (const_vals.Ln)
    {
        this->e_ny = integer(0);
        this->num_coeffs = this->num_coeffs - 1;
    }
    if (const_vals.Rn)
    {
        this->f_ny = integer(0);
        this->num_coeffs = this->num_coeffs - 1;
    }
    if (const_vals.Tn)
    {
        this->d_hx = integer(0);
        this->num_coeffs = this->num_coeffs - 1;
    }
    if (const_vals.Bn)
    {
        this->c_hx = integer(0);
        this->num_coeffs = this->num_coeffs - 1;
    }
}
