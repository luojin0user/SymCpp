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

BasicCase::BasicCase(int idx,
                     const RCP<const Symbol> &xl,
                     const RCP<const Symbol> &xr,
                     const RCP<const Symbol> &yl,
                     const RCP<const Symbol> &yt,
                     int H_max, int N_max,
                     const RCP<const Symbol> &mu_r)
{
    this->idx = idx;

    // syms x y
    this->x = symbol("x");
    this->y = symbol("y");

    std::string suffix = "_" + std::to_string(idx);

    // --- 公共编号变量 ---
    // MATLAB 逻辑: 优先使用传入的 xl, xr...
    this->xr = xr;
    this->xl = xl;
    this->yt = yt;
    this->yl = yl;

    // obj.tau_x = obj.xr - obj.xl;
    this->tau_x = sub(this->xr, this->xl);
    this->tau_y = sub(this->yt, this->yl);

    // --- 公共求和指标 ---
    this->h = symbol("h" + suffix);
    this->n = symbol("n" + suffix);

    this->H_max = H_max;
    this->N_max = N_max;

    // --- 公共参数 ---
    // obj.beta_h = obj.h * pi / obj.tau_x;
    this->beta_h = div(mul(this->h, pi), this->tau_x);

    // obj.lambda_n = obj.n * pi / obj.tau_y;
    this->lambda_n = div(mul(this->n, pi), this->tau_y);

    // obj.mu_0 = 4*pi*1e-7;
    this->mu_0 = mul(mul(integer(4), pi), real_double(1e-7));
    this->mu_r = mu_r;

    // --- 边界积分项 ---
    this->c_0x = symbol("c_0x" + suffix);
    this->d_0x = symbol("d_0x" + suffix);
    this->c_hx = symbol("c_hx" + suffix);
    this->d_hx = symbol("d_hx" + suffix);
    this->e_ny = symbol("e_ny" + suffix);
    this->f_ny = symbol("f_ny" + suffix);

    // 初始化 num_coeffs (MATLAB中未在构造函数显式初始化的部分属性，需结合后续逻辑)
    this->num_coeffs = 6;
}

void BasicCase::apply_boundaries(bool Ln, bool Rn, bool Tn, bool Bn)
{
    this->Ln = Ln;
    this->Rn = Rn;
    this->Tn = Tn;
    this->Bn = Bn;

    if (Ln)
    {
        this->e_ny = integer(0);
        this->num_coeffs = this->num_coeffs - 1;
    }
    if (Rn)
    {
        this->f_ny = integer(0);
        this->num_coeffs = this->num_coeffs - 1;
    }
    if (Tn)
    {
        this->d_hx = integer(0);
        this->num_coeffs = this->num_coeffs - 1;
    }
    if (Bn)
    {
        this->c_hx = integer(0);
        this->num_coeffs = this->num_coeffs - 1;
    }
}

void BasicCase::gen_solution_func(bool Ln, bool Rn, bool Tn, bool Bn)
{
    throw std::runtime_error("子类必须实现 gen_solution_func 方法");
}

void BasicCase::gen_coefficient_func()
{
    throw std::runtime_error("子类必须实现 gen_coefficient_func 方法");
}
