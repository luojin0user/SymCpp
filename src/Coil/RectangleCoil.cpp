#include "RectangleCoil.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

// SymEngine 具体实现头文件
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/integer.h>
#include <symengine/real_double.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/functions.h>
#include <symengine/constants.h>
#include <symengine/eval_double.h>

// 在 cpp 中使用 using namespace
using SymEngine::acos;
using SymEngine::add;
using SymEngine::Basic;
using SymEngine::cos;
using SymEngine::div;
using SymEngine::integer;
using SymEngine::mul;
using SymEngine::pi;
using SymEngine::pow;
using SymEngine::RCP;
using SymEngine::real_double;
using SymEngine::sin;
using SymEngine::sqrt;
using SymEngine::sub;
using SymEngine::symbol;
using SymEngine::Symbol;

RectangleCoil::RectangleCoil(double outer_len_x, double outer_len_y, int Nt_per_layer, int layer, double wire_diameter, double Ir)
{
    this->outer_len_x = outer_len_x;
    this->outer_len_y = outer_len_y;
    this->Nt_per_layer = Nt_per_layer;
    this->layer = layer;
    this->Ir = Ir;
    this->wire_diameter = wire_diameter;

    // obj.inner_len_x = outer_len_x - Nt_per_layer * wire_diameter * 2;
    this->inner_len_x = outer_len_x - Nt_per_layer * wire_diameter * 2.0;
    this->inner_len_y = outer_len_y - Nt_per_layer * wire_diameter * 2.0;

    // obj.Nt = Nt_per_layer * layer;
    this->Nt = Nt_per_layer * layer;
    // obj.thick = wire_diameter * layer;
    this->thick = wire_diameter * layer;

    // 初始化部分默认值
    this->x_loc = 0;
    this->y_loc = 0;
    this->z_loc = 0;
    this->mu0 = 4.0 * 3.141592653589793 * 1e-7;
}

void RectangleCoil::gen_all_regions()
{
    // plain = AllRegions(); ...
    // 因缺少依赖，保留空实现或打印警告
    std::cout << "[RectangleCoil::gen_all_regions] Warning: AllRegions class missing. Function skipped." << std::endl;
}

void RectangleCoil::set_Rcoil_loc(double x, double y, double z)
{
    this->x_loc = x;
    this->y_loc = y;
    this->z_loc = z;

    this->x_cl1 = x - this->outer_len_x / 2.0;
    this->x_cl2 = x - this->inner_len_x / 2.0;
    this->x_cr1 = x + this->inner_len_x / 2.0;
    this->x_cr2 = x + this->outer_len_x / 2.0;

    this->y_cl1 = y - this->outer_len_y / 2.0;
    this->y_cl2 = y - this->inner_len_y / 2.0;
    this->y_cr1 = y + this->inner_len_y / 2.0;
    this->y_cr2 = y + this->outer_len_y / 2.0;

    this->m_xl = (this->x_cl1 + this->x_cl2) / 2.0;
    this->m_xr = (this->x_cr1 + this->x_cr2) / 2.0;
    this->m_yl = (this->y_cl1 + this->y_cl2) / 2.0;
    this->m_yr = (this->y_cr1 + this->y_cr2) / 2.0;

    this->z_b = z;
    this->z_t = z + this->thick;

    // obj.Jr_xl = obj.Nt * obj.Ir / (obj.thick .* (obj.outer_len_x - obj.inner_len_x) ./ 2);
    double area_x = this->thick * (this->outer_len_x - this->inner_len_x) / 2.0;
    this->Jr_xl = (this->Nt * this->Ir) / area_x;
    this->Jr_xr = -this->Jr_xl;

    // obj.Jr_yl = obj.Nt * obj.Ir / (obj.thick .* (obj.outer_len_y - obj.inner_len_y) ./ 2);
    double area_y = this->thick * (this->outer_len_y - this->inner_len_y) / 2.0;
    this->Jr_yl = (this->Nt * this->Ir) / area_y;
    this->Jr_yr = -this->Jr_yl;
}

void RectangleCoil::set_calculate_area(double xl, double xr, double yl, double yr, double zb, double zt, double mu0)
{
    this->cal_xl = xl;
    this->cal_xr = xr;
    this->cal_yl = yl;
    this->cal_yr = yr;
    this->cal_zb = zb;
    this->cal_zt = zt;

    this->mu0 = mu0;

    this->factor_x();
    this->factor_y();
}

void RectangleCoil::factor_x()
{
    RCP<const Symbol> x = symbol("x");
    RCP<const Symbol> y = symbol("y");
    RCP<const Symbol> z = symbol("z");

    // 将 double 参数转换为 SymEngine 数值对象
    RCP<const Basic> mu_0_s = real_double(this->mu0);
    RCP<const Basic> Ip_s = real_double(this->Ir);
    RCP<const Basic> Zc_s = real_double(this->thick);
    RCP<const Basic> m_xl_s = real_double(this->m_xl);
    RCP<const Basic> m_xr_s = real_double(this->m_xr);
    RCP<const Basic> y_cl1_s = real_double(this->y_cl1);
    RCP<const Basic> y_cr2_s = real_double(this->y_cr2);

    // g1 = sqrt((x - obj.m_xl)^2 + (z - Zc)^2);
    RCP<const Basic> g1 = sqrt(add(pow(sub(x, m_xl_s), integer(2)), pow(sub(z, Zc_s), integer(2))));

    // g2 = sqrt((x - obj.m_xr)^2 + (z - Zc)^2);
    RCP<const Basic> g2 = sqrt(add(pow(sub(x, m_xr_s), integer(2)), pow(sub(z, Zc_s), integer(2))));

    // g3 = obj.y_cl1 - y;
    RCP<const Basic> g3 = sub(y_cl1_s, y);

    // g4 = obj.y_cr2 - y;
    RCP<const Basic> g4 = sub(y_cr2_s, y);

    RCP<const Basic> two_pi = mul(integer(2), pi);
    RCP<const Basic> g1_2pi = mul(g1, two_pi);
    RCP<const Basic> g2_2pi = mul(g2, two_pi);

    // alpha1 = acos((x - obj.m_xl) ./ g1);
    RCP<const Basic> alpha1 = acos(div(sub(x, m_xl_s), g1));

    // alpha2 = acos((obj.m_xr - x) ./ g2);
    RCP<const Basic> alpha2 = acos(div(sub(m_xr_s, x), g2));

    // B_inf_x 计算
    RCP<const Basic> term1 = div(mul(mu_0_s, Ip_s), g1_2pi);
    RCP<const Basic> term2 = div(mul(mu_0_s, Ip_s), g2_2pi);

    // cos(pi - alpha1 - alpha2)
    RCP<const Basic> cos_term = cos(sub(sub(pi, alpha1), alpha2));

    // B_inf_x = sqrt(term1^2 + term2^2 - 2*term1*term2*cos(...))
    RCP<const Basic> B_inf_x = sqrt(sub(add(pow(term1, integer(2)), pow(term2, integer(2))),
                                        mul(mul(integer(2), mul(term1, term2)), cos_term)));

    // B_f_1
    // coeff = mu_0 * Ip / (2 * g1_2pi) -> 注意 MATLAB 代码中是 mu_0 .* Ip ./ (2 .* g1_2pi) ?
    // 检查 MATLAB: mu_0 .* Ip ./ g1_2pi 是用于 B_inf_x 的.
    // 对于 B_f_1: mu_0 .* Ip ./ (2 .* g1_2pi) -> 这里的 g1_2pi 已经是 2*pi*g1 了
    // 看起来 MATLAB 里的命名 g1_2pi 是 2*pi*g1.
    // B_f_1 = ... ./ (2 .* g1_2pi) 意味着除以 4*pi*g1 ?
    // 让我们照搬 MATLAB 公式： mu_0 .* Ip ./ (2 .* g1_2pi)
    RCP<const Basic> bf1_coeff = div(mul(mu_0_s, Ip_s), mul(integer(2), g1_2pi));
    RCP<const Basic> bf1_p1 = div(g4, sqrt(add(pow(g1, integer(2)), pow(g4, integer(2)))));
    RCP<const Basic> bf1_p2 = div(g3, sqrt(add(pow(g1, integer(2)), pow(g3, integer(2)))));
    RCP<const Basic> B_f_1 = mul(bf1_coeff, sub(bf1_p1, bf1_p2));

    // B_f_2
    RCP<const Basic> bf2_coeff = div(mul(mu_0_s, Ip_s), mul(integer(2), g2_2pi));
    RCP<const Basic> bf2_p1 = div(g4, sqrt(add(pow(g2, integer(2)), pow(g4, integer(2)))));
    RCP<const Basic> bf2_p2 = div(g3, sqrt(add(pow(g2, integer(2)), pow(g3, integer(2)))));
    RCP<const Basic> B_f_2 = mul(bf2_coeff, sub(bf2_p1, bf2_p2));

    // B_f_x
    RCP<const Basic> B_f_x = sqrt(sub(add(pow(B_f_1, integer(2)), pow(B_f_2, integer(2))),
                                      mul(mul(integer(2), mul(B_f_1, B_f_2)), cos_term)));

    // obj.fx = matlabFunction(expr...)
    this->fx = div(B_f_x, B_inf_x);
}

double RectangleCoil::factor_x_real(double x_val, double y_val, double z_val)
{
    if (this->fx.is_null())
    {
        this->factor_x();
    }

    SymEngine::map_basic_basic subs_map;
    subs_map[symbol("x")] = real_double(x_val);
    subs_map[symbol("y")] = real_double(y_val);
    subs_map[symbol("z")] = real_double(z_val);

    auto res = this->fx->subs(subs_map);
    return SymEngine::eval_double(*res);
}

void RectangleCoil::factor_y()
{
    RCP<const Symbol> x = symbol("x");
    RCP<const Symbol> y = symbol("y");
    RCP<const Symbol> z = symbol("z");

    RCP<const Basic> mu_0_s = real_double(this->mu0);
    RCP<const Basic> Ip_s = real_double(this->Ir);
    RCP<const Basic> Zc_s = real_double(this->thick);
    RCP<const Basic> m_yl_s = real_double(this->m_yl);
    RCP<const Basic> m_yr_s = real_double(this->m_yr);
    RCP<const Basic> x_cl1_s = real_double(this->x_cl1);
    RCP<const Basic> x_cr2_s = real_double(this->x_cr2);

    // 注意：与 factor_x 相比，x 和 y 的角色互换，基于 MATLAB 源码逻辑
    // g1 = sqrt((y - obj.m_yl)^2 + (z - Zc)^2);
    RCP<const Basic> g1 = sqrt(add(pow(sub(y, m_yl_s), integer(2)), pow(sub(z, Zc_s), integer(2))));

    // g2 = sqrt((y - obj.m_yr)^2 + (z - Zc)^2);
    RCP<const Basic> g2 = sqrt(add(pow(sub(y, m_yr_s), integer(2)), pow(sub(z, Zc_s), integer(2))));

    // g3 = obj.x_cl1 - x;
    RCP<const Basic> g3 = sub(x_cl1_s, x);

    // g4 = obj.x_cr2 - x;
    RCP<const Basic> g4 = sub(x_cr2_s, x);

    RCP<const Basic> two_pi = mul(integer(2), pi);
    RCP<const Basic> g1_2pi = mul(g1, two_pi);
    RCP<const Basic> g2_2pi = mul(g2, two_pi);

    // alpha1 = acos((y - obj.m_yl) ./ g1);
    RCP<const Basic> alpha1 = acos(div(sub(y, m_yl_s), g1));

    // alpha2 = acos((obj.m_yr - y) ./ g2);
    RCP<const Basic> alpha2 = acos(div(sub(m_yr_s, y), g2));

    // B_inf_y calculation
    RCP<const Basic> term1 = div(mul(mu_0_s, Ip_s), g1_2pi);
    RCP<const Basic> term2 = div(mul(mu_0_s, Ip_s), g2_2pi);
    RCP<const Basic> cos_term = cos(sub(sub(pi, alpha1), alpha2));

    RCP<const Basic> B_inf_y = sqrt(sub(add(pow(term1, integer(2)), pow(term2, integer(2))),
                                        mul(mul(integer(2), mul(term1, term2)), cos_term)));

    // B_f_1
    RCP<const Basic> bf1_coeff = div(mul(mu_0_s, Ip_s), mul(integer(2), g1_2pi));
    RCP<const Basic> bf1_p1 = div(g4, sqrt(add(pow(g1, integer(2)), pow(g4, integer(2)))));
    RCP<const Basic> bf1_p2 = div(g3, sqrt(add(pow(g1, integer(2)), pow(g3, integer(2)))));
    RCP<const Basic> B_f_1 = mul(bf1_coeff, sub(bf1_p1, bf1_p2));

    // B_f_2
    RCP<const Basic> bf2_coeff = div(mul(mu_0_s, Ip_s), mul(integer(2), g2_2pi));
    RCP<const Basic> bf2_p1 = div(g4, sqrt(add(pow(g2, integer(2)), pow(g4, integer(2)))));
    RCP<const Basic> bf2_p2 = div(g3, sqrt(add(pow(g2, integer(2)), pow(g3, integer(2)))));
    RCP<const Basic> B_f_2 = mul(bf2_coeff, sub(bf2_p1, bf2_p2));

    // B_f_y
    RCP<const Basic> B_f_y = sqrt(sub(add(pow(B_f_1, integer(2)), pow(B_f_2, integer(2))),
                                      mul(mul(integer(2), mul(B_f_1, B_f_2)), cos_term)));

    this->fy = div(B_f_y, B_inf_y);
}