#include "Boundarys.hpp"
#include "Region.hpp"
#include "BasicCase.hpp"
#include "Case1.hpp" // For dynamic_cast or member access
#include "Case2.hpp"

#include <iostream>

// SymEngine includes
#include <symengine/functions.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/subs.h>

#include "EnumTypes.h"

using SymEngine::Basic;
using SymEngine::integer;
using SymEngine::map_basic_basic;
using SymEngine::RCP;
using SymEngine::real_double;
using SymEngine::Symbol;
using SymEngine::symbol;

Boundarys::Boundarys(const Region *impl, int all_region_num)
    : region(impl), region_impl(impl->impl) // 初始化引用
{
    this->all_region_num = all_region_num;
}

void Boundarys::cal_BC(bool Ln, bool Rn, bool Tn, bool Bn)
{
    int BCfuncs_loc_num = 1;

    switch (this->bc_type)
    {
    case BC_TYPE::BBAA:
        // 这里需要注意可能出现的分段函数的形式，如果出现分段函数，一定是B连续，所以需要
        // 目前而言，出现分段函数只会在上下边界为BB的情况出现
        // 当出现分段时，采用的储存方法为BC{i}{j}，BC{1}{j}是上边的第一个分段区域，BC{2}{j}是上边的第二个分段区域，以此类推
        // 这样表示的方法和论文中的类似
        // 上边界
        process_boundary_BB(!Tn, this->top, region_impl->y, region_impl->yt, T_funcs, T_ESfuncs);
        // 下边界
        process_boundary_BB(!Bn, this->bottom, region_impl->y, region_impl->yl, B_funcs, B_ESfuncs);

        // 左右边界由于是A，需要乘以region_impl.lambdaN
        // 如果上下边界为A，则需要乘以region_impl.betaH
        process_boundary_AA(!Ln, this->left, region_impl->x, region_impl->xl, 2, L_funcs, L_ESfuncs);
        process_boundary_AA(!Rn, this->right, region_impl->x, region_impl->xr, 2, R_funcs, R_ESfuncs);
        break;

    case BC_TYPE::AAAA:
        process_boundary_AA(!Tn, this->top, region_impl->y, region_impl->yt, 1, T_funcs, T_ESfuncs);
        process_boundary_AA(!Bn, this->bottom, region_impl->y, region_impl->yl, 1, B_funcs, B_ESfuncs);

        process_boundary_AA(!Ln, this->left, region_impl->x, region_impl->xl, 2, L_funcs, L_ESfuncs);
        process_boundary_AA(!Rn, this->right, region_impl->x, region_impl->xr, 2, R_funcs, R_ESfuncs);
        break;

    case BC_TYPE::AABB:
        process_boundary_AA(!Tn, this->top, region_impl->y, region_impl->yt, 1, T_funcs, T_ESfuncs);
        process_boundary_AA(!Bn, this->bottom, region_impl->y, region_impl->yl, 1, B_funcs, B_ESfuncs);

        process_boundary_BB(!Ln, this->left, region_impl->x, region_impl->xl, L_funcs, L_ESfuncs);
        process_boundary_BB(!Rn, this->right, region_impl->x, region_impl->xr, R_funcs, R_ESfuncs);
        break;
    }
}

void Boundarys::cal_ES(bool Ln, bool Rn, bool Tn, bool Bn)
{
    // Empty in MATLAB source provided
}

Boundary_Funcs Boundarys::genAA(int right_idx)
{
    // 首先需要找到对应的方程，如果对应的方程中间包含多个c d e f，需要一一进行判断处理，然后送入数组中，和这个区域的边界情况的方程匹配
    // 如果由多个，送入数组中，由region类对这个数组进行处理
    Boundary_Funcs res;

    const auto &edge_region = *this->region->all_regions->at(right_idx); // 获取邻接区域对象
    const auto &edge_impl = *edge_region.impl;
    bool has_cd0x = (edge_region.case_type == CaseType::FerriteCurrent);

    RCP<const Basic> F1;

    F1 = edge_impl.A_zx_expr;
    if (has_cd0x)
    {
        // 首先处理c_0x
        map_basic_basic d;
        d[edge_impl.c_hx] = integer(0);
        d[edge_impl.d_hx] = integer(0);
        d[edge_impl.c_0x] = integer(1);
        d[edge_impl.d_0x] = integer(0);
        res.funcs[0] = F1->subs(d);

        d[edge_impl.c_0x] = integer(0);
        d[edge_impl.d_0x] = integer(1);
        res.funcs[2] = F1->subs(d);
    }

    // --- 2 & 4: c_hx / d_hx ---
    if (!edge_impl.const_vals.Bn)
    {
        map_basic_basic d;
        d[edge_impl.c_hx] = integer(1);
        d[edge_impl.d_hx] = integer(0);
        if (has_cd0x)
        {
            d[edge_impl.c_0x] = integer(0);
            d[edge_impl.d_0x] = integer(0);
        }

        res.funcs[1] = F1->subs(d);
    }
    if (!edge_impl.const_vals.Tn)
    {
        map_basic_basic d;
        d[edge_impl.c_hx] = integer(0);
        d[edge_impl.d_hx] = integer(1);
        if (has_cd0x)
        {
            d[edge_impl.c_0x] = integer(0);
            d[edge_impl.d_0x] = integer(0);
        }

        res.funcs[3] = F1->subs(d);
    }

    // --- 5 & 6: e_ny / f_ny ---
    F1 = edge_impl.A_zy_expr;
    if (!edge_impl.const_vals.Ln)
    {
        map_basic_basic d;
        d[edge_impl.e_ny] = integer(1);
        d[edge_impl.f_ny] = integer(0);
        res.funcs[4] = F1->subs(d);
    }
    if (!edge_impl.const_vals.Rn)
    {
        map_basic_basic d;
        d[edge_impl.e_ny] = integer(0);
        d[edge_impl.f_ny] = integer(1);
        res.funcs[5] = F1->subs(d);
    }

    res.idx = edge_impl.const_vals.idx;
    res.rect = edge_impl.const_vals.rect;

    return res;
}

Boundary_Funcs Boundarys::genBB(int right_idx, int cd_or_ef)
{
    // 这里记得乘以2个的mu_0的系数
    // 首先需要找到对应的方程，如果对应的方程中间包含多个c d e f，需要一一进行判断处理，然后送入数组中，和这个区域的边界情况的方程匹配
    // 如果由多个，送入数组中，由region类对这个数组进行处理
    Boundary_Funcs res;

    // 邻接区域的 A_z，
    // 当是邻接区域的，直接把这个邻接区域的方程A_z送给对应的top
    // 除此之外，还要令对应的方程的变量值为边界值
    // 对于一个已经确定的区域，其e、f参数由B_y决定，其中某一个(e)的参数c d由B_y_x决定，e f由B_y_y决定
    const auto &edge_region = *this->region->all_regions->at(right_idx); // 获取邻接区域对象
    const auto &edge_impl = *edge_region.impl;
    bool has_cd0x = (edge_region.case_type == CaseType::FerriteCurrent);

    RCP<const Basic> mu_ratio = SymEngine::div(this->region_impl->mu_r, edge_impl.mu_r); // 乘以 mu_r / mu_r

    RCP<const Basic> F1;

    // c_0x / d_0x
    if (has_cd0x)
    {
        F1 = this->getBB_func(edge_impl, 1, cd_or_ef); // 1 = x direction for linear term check?
        // 首先处理c_0x
        map_basic_basic d;
        d[edge_impl.c_hx] = integer(0);
        d[edge_impl.d_hx] = integer(0);
        d[edge_impl.c_0x] = integer(1);
        d[edge_impl.d_0x] = integer(0);
        res.funcs[0] = F1->subs(d);

        d[edge_impl.c_0x] = integer(0);
        d[edge_impl.d_0x] = integer(1);
        res.funcs[2] = F1->subs(d);

        res.funcs[0] = SymEngine::mul(res.funcs[0], mu_ratio);
        res.funcs[2] = SymEngine::mul(res.funcs[2], mu_ratio);
    }

    // c_hx
    if (!edge_impl.const_vals.Bn)
    {
        F1 = this->getBB_func(edge_impl, 1, cd_or_ef);
        map_basic_basic d;
        d[edge_impl.c_hx] = integer(1);
        d[edge_impl.d_hx] = integer(0);
        if (has_cd0x)
        {
            d[edge_impl.c_0x] = integer(0);
            d[edge_impl.d_0x] = integer(0);
        }

        F1 = F1->subs(d);

        res.funcs[1] = SymEngine::mul(F1, mu_ratio);
    }
    // d_hx
    if (!edge_impl.const_vals.Tn)
    {
        F1 = this->getBB_func(edge_impl, 1, cd_or_ef);
        map_basic_basic d;
        d[edge_impl.c_hx] = integer(0);
        d[edge_impl.d_hx] = integer(1);
        if (has_cd0x)
        {
            d[edge_impl.c_0x] = integer(0);
            d[edge_impl.d_0x] = integer(0);
        }
        res.funcs[3] = F1->subs(d);

        res.funcs[3] = SymEngine::mul(res.funcs[3], mu_ratio);
    }

    // e_ny / f_ny
    if (!edge_impl.const_vals.Ln)
    {
        F1 = this->getBB_func(edge_impl, 2, cd_or_ef);
        map_basic_basic d;
        d[edge_impl.e_ny] = integer(1);
        d[edge_impl.f_ny] = integer(0);
        res.funcs[4] = F1->subs(d);

        res.funcs[4] = SymEngine::mul(res.funcs[4], mu_ratio);
    }
    if (!edge_impl.const_vals.Rn)
    {
        F1 = this->getBB_func(edge_impl, 2, cd_or_ef);
        map_basic_basic d;
        d[edge_impl.e_ny] = integer(0);
        d[edge_impl.f_ny] = integer(1);
        res.funcs[5] = F1->subs(d);

        res.funcs[5] = SymEngine::mul(res.funcs[5], mu_ratio);
    }

    res.idx = edge_impl.const_vals.idx;
    res.rect = edge_impl.const_vals.rect;

    return res;
}

RCP<const Basic> Boundarys::getBB_func(const BasicCase &edge_impl, int x_or_y, int cd_or_ef)
{
    RCP<const Basic> F1;
    if (x_or_y == 1) // 上下边界为x，左右为y
    {                // x
        if (cd_or_ef == 1)
            F1 = edge_impl.B_x_x; // 只含c d项
        else
            F1 = edge_impl.B_y_x;
    }
    else
    { // y
        if (cd_or_ef == 1)
            F1 = edge_impl.B_x_y;
        else
            F1 = edge_impl.B_y_y;
    }
    return F1;
}

void Boundarys::process_boundary_BB(
    bool BC_exist,                            // 边界是否存在 !Bn，如果为true则存在边界
    const std::vector<int> &boundary_indices, // 当前处理的的边界
    const SymEngine::RCP<const SymEngine::Basic> &sub_src,
    const SymEngine::RCP<const SymEngine::Basic> &sub_dst,
    std::vector<Boundary_Funcs> &funcs,                            // 需要返回的边界方程
    std::vector<SymEngine::RCP<const SymEngine::Basic>> &ES_funcs) // 可能存在的ES方程
{

    auto subs_func_with_map = [](Boundary_Funcs &funcs, const RCP<const Basic> &val, const RCP<const Basic> &x)
    {
        map_basic_basic d;
        d[val] = x;

        for (auto &i : funcs.funcs)
            if (!i.is_null())
                i = i->subs(d);
    };

    auto handle_ES = [&](const BasicCase &edge_impl, const RCP<const Basic> &val, const RCP<const Basic> &x)
    {
        auto ES = ((Case2 &)(edge_impl)).B_x_P;
        map_basic_basic d;
        d[val] = x;
        ES = ES->subs(d);

        return ES;
    };

    if (!BC_exist || boundary_indices.empty())
    {
        funcs.clear();
        return;
    }

    for (int idx : boundary_indices)
    {
        // 生成边界方程
        Boundary_Funcs eq = this->genBB(idx, 1);
        subs_func_with_map(eq, sub_src, sub_dst);
        funcs.push_back(eq);

        // 处理 ES
        if (this->region->impl->ES_regions[idx])
        {
            const auto &edge_impl = *this->region->all_regions->at(idx)->impl;
            ES_funcs.push_back(handle_ES(edge_impl, sub_src, sub_dst));
        }
    }
}

void Boundarys::process_boundary_AA(
    bool BC_exist,                            // 边界是否存在 !Bn，如果为true则存在边界
    const std::vector<int> &boundary_indices, // 当前处理的的边界
    const SymEngine::RCP<const SymEngine::Basic> &sub_src,
    const SymEngine::RCP<const SymEngine::Basic> &sub_dst,
    int tb_or_lr,                                                  // 当前处理的是上下(tb,1)或者左右(lr,2)
    std::vector<Boundary_Funcs> &funcs,                            // 需要返回的边界方程
    std::vector<SymEngine::RCP<const SymEngine::Basic>> &ES_funcs) // 可能存在的ES方程
{

    auto subs_func_with_map = [&](Boundary_Funcs &funcs, const RCP<const Basic> &val, const RCP<const Basic> &x)
    {
        map_basic_basic d;
        d[val] = x;

        SymEngine::RCP<const SymEngine::Basic> xm;
        if (tb_or_lr == 1)
            xm = region_impl->beta_h;
        else
            xm = region_impl->lambda_n;

        for (auto &i : funcs.funcs)
            if (!i.is_null())
            {
                i = i->subs(d);
                i = SymEngine::mul(i, xm);
            }
    };

    auto handle_ES = [&](const BasicCase &edge_impl, const RCP<const Basic> &val, const RCP<const Basic> &x)
    {
        SymEngine::RCP<const SymEngine::Basic> xm;
        if (tb_or_lr == 1)
            xm = region_impl->beta_h;
        else
            xm = region_impl->lambda_n;

        auto ES = ((Case2 &)(edge_impl)).B_x_P;
        map_basic_basic d;
        d[val] = x;
        ES = ES->subs(d);

        ES = SymEngine::mul(ES, xm);
        return ES;
    };

    if (!BC_exist || boundary_indices.empty())
    {
        funcs.clear();
        return;
    }

    for (int idx : boundary_indices)
    {
        // 生成边界方程
        Boundary_Funcs eq = this->genAA(idx);
        subs_func_with_map(eq, sub_src, sub_dst);
        funcs.push_back(eq);

        // 处理 ES
        if (this->region->impl->ES_regions[idx])
        {
            const auto &edge_impl = *this->region->all_regions->at(idx)->impl;
            ES_funcs.push_back(handle_ES(edge_impl, sub_src, sub_dst));
        }
    }
}

void Boundarys::rtn_BC_funcs(std::vector<Boundary_Funcs> &T_funcs,
                             std::vector<Boundary_Funcs> &B_funcs,
                             std::vector<Boundary_Funcs> &L_funcs,
                             std::vector<Boundary_Funcs> &R_funcs)
{
    T_funcs.assign(this->T_funcs.begin(), this->T_funcs.end());
    B_funcs.assign(this->B_funcs.begin(), this->B_funcs.end());
    L_funcs.assign(this->L_funcs.begin(), this->L_funcs.end());
    R_funcs.assign(this->R_funcs.begin(), this->R_funcs.end());
}

void Boundarys::rtn_ES_funcs(std::vector<SymEngine::RCP<const SymEngine::Basic>> &T_ESfuncs,
                             std::vector<SymEngine::RCP<const SymEngine::Basic>> &B_ESfuncs,
                             std::vector<SymEngine::RCP<const SymEngine::Basic>> &L_ESfuncs,
                             std::vector<SymEngine::RCP<const SymEngine::Basic>> &R_ESfuncs)
{
    T_ESfuncs.assign(this->T_ESfuncs.begin(), this->T_ESfuncs.end());
    B_ESfuncs.assign(this->B_ESfuncs.begin(), this->B_ESfuncs.end());
    L_ESfuncs.assign(this->L_ESfuncs.begin(), this->L_ESfuncs.end());
    R_ESfuncs.assign(this->R_ESfuncs.begin(), this->R_ESfuncs.end());
}