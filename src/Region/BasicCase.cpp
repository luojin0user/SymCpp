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

void BasicCase::gen_integral(const std::vector<Boundary_Funcs> &BC_func,
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