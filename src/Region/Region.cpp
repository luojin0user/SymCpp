#include "Region.hpp"
#include "RegionsInput.hpp"
#include <iostream>
#include <unordered_map>
#include "EnumTypes.h"

#include "AlleyAir.hpp"
#include "BTAir.hpp"
#include "FerriteCurrent.hpp"
#include "NormalAir.hpp"

// 使用命名空间以便于访问 SymEngine 类
using SymEngine::Basic;
using SymEngine::RCP;
using SymEngine::Symbol;

Region::Region(int idx, CaseType type,
               const Rect &area,
               BC_TYPE bc_type,

               const std::vector<int> &top,
               const std::vector<int> &bottom,
               const std::vector<int> &left,
               const std::vector<int> &right,

               const std::vector<int> &ES_regions_idx,

               int H_max, int N_max,
               float mu_r,
               float J_r,
               int all_regions_num)
{
    // L/R/T/B 为 bool 值，表示边界是否存在，如果为True，说明这个边界不存在(即有邻居)
    // obj.Ln = isempty(left);
    Region_Consts c;
    c.idx = idx;
    c.rect = area;
    c.H_max = H_max;
    c.N_max = N_max;
    c.mu_r = mu_r;
    c.J_r = J_r;
    c.all_regions_num = all_regions_num;

    c.set_Boundarys(top.empty(), bottom.empty(), left.empty(), right.empty());

    this->case_type = type;

    // Factory Pattern Implementation
    switch (type)
    {
    case CaseType::AlleyAir:
        this->impl = std::make_unique<AlleyAir>(c);
        break;
    case CaseType::BTAir: // NormalAir logic in MATLAB comments?
        this->impl = std::make_unique<BTAir>(c);
        break;
    case CaseType::NormalAir:
        this->impl = std::make_unique<NormalAir>(c);
        break;
    case CaseType::FerriteCurrent:
        this->impl = std::make_unique<FerriteCurrent>(c);
        break;
    case CaseType::Aluminum:
        // this->impl = std::make_unique<Aluminum>(idx, this->Ln, this->Rn, this->Tn, this->Bn, H_max, N_max, mu_r);
        break;
    default:
        throw std::runtime_error("Unsupported case type");
    }

    this->impl->ES_regions.assign(all_regions_num, false);

    for (int es_idx : ES_regions_idx)
    {
        if (es_idx >= 0 && es_idx < all_regions_num)
        {
            this->impl->ES_regions[es_idx] = true;
        }
    }

    this->all_regions_num = all_regions_num;

    this->boundarys = std::make_unique<Boundarys>(this, this->all_regions_num);
    this->set_boundary(top, bottom, left, right, bc_type);
}

void Region::set_boundary(const std::vector<int> &top,
                          const std::vector<int> &bottom,
                          const std::vector<int> &left,
                          const std::vector<int> &right,
                          BC_TYPE bc_type)
{
    this->boundarys->bottom = bottom;
    this->boundarys->top = top;
    this->boundarys->left = left;
    this->boundarys->right = right;
    this->boundarys->bc_type = bc_type;

    // obj.all_edge_regions = false(obj.all_regions_num, 1);
    all_edge_regions.resize(all_regions_num);
    for (int i = 0; i < this->all_regions_num; i++)
        this->all_edge_regions[i] = {false, -1};

    // obj.all_edge_regions(regions) = true;
    for (int r_idx : bottom)
        this->all_edge_regions[r_idx] = {true, 1};
    for (int r_idx : top)
        this->all_edge_regions[r_idx] = {true, 3};
    for (int r_idx : left)
        this->all_edge_regions[r_idx] = {true, 4};
    for (int r_idx : right)
        this->all_edge_regions[r_idx] = {true, 5};

    cal_BCxx_loc();
}

// 计算这个区域内部的方程，这里直接返回这个方程，计算过程由set_boundary中的obj.impl.gen_region完成
void Region::get_region_solution_func()
{
    // obj.impl.gen_solution_func();
    // Using the overloaded version or default based on previous impl
    this->impl->gen_solution_func();

    // region_func = {obj.impl.eq_A_z, obj.impl.eq_B_x};
}

// 计算这个区域系数的方程，这里是在所有区域内部方程计算完成之后调用的
void Region::gen_region_coefficient_func(const std::vector<std::unique_ptr<Region>> *all_regions)
{
    // 首先需要计算这个区域的边界的情况，需要使用boundays.cal_BC方法，边界方程的计算需要用到所有的区域
    this->all_regions = all_regions;

    this->boundarys->cal_BC(this->impl->const_vals.Ln, this->impl->const_vals.Rn,
                            this->impl->const_vals.Tn, this->impl->const_vals.Bn);

    // 这里会返回所有边界方程的数据
    this->boundarys->rtn_BC_funcs(this->impl->T_funcs, this->impl->B_funcs, this->impl->L_funcs, this->impl->R_funcs);
    this->boundarys->rtn_ES_funcs(this->impl->T_ESfuncs, this->impl->B_ESfuncs, this->impl->L_ESfuncs, this->impl->R_ESfuncs);

    this->impl->gen_coefficient_func();

    // 这里需要组合所有c0 c d0 d e f方程，方程的定义在BasicCase中
    // std::array<std::vector<std::array<Integral_Func, 6>>, 6> eq_BC;
    // 每个Integral_Func结构体里面储存了邻接区域的编号信息
    // 直接被调用
    return;
}

void Region::cal_BCxx_loc()
{
    int curr_loc = -1;
    int curr_region = 0;

    // 用于指示边界位置（BTLR 0123）对应在BCxx中的位置
    // 理论上的位置就是c0 c d0 d e f，如果哪个没有出现那么就删去
    // 是否出现的flag存在coeffs_exists中
    std::unordered_map<int, int> row_exists_idx;
    int counter = 0;
    for (int j = 0; j < this->impl->coeffs_exists.size(); j++)
    {
        if (this->impl->coeffs_exists[j])
        {
            // 如果这个边界存在，记录下这个的位置
            row_exists_idx[j] = counter++;
        }
    }

    // 对于所有的边界
    for (int i = 0; i < all_edge_regions.size(); i++)
    {
        // 存在这个边界，如果没有删去放在 .second行
        if (all_edge_regions[i].first)
        {
            auto find_row_idx = row_exists_idx.find(all_edge_regions[i].second);
            if (find_row_idx != row_exists_idx.end())
            {
                // 如果这个边界存在，记录下这个的位置
                BCxx_loc[all_edge_regions[i].second] = find_row_idx->second; // 顺序
            }

            // 处理有源项目
            if (this->case_type == CaseType::FerriteCurrent &&
                (all_edge_regions[i].second == 1 || all_edge_regions[i].second == 3))
            {
                BCxx_loc[all_edge_regions[i].second - 1] = find_row_idx->second - 1;
            }
        }
    }
}