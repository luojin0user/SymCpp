#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

// SymEngine
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/real_double.h>
#include <symengine/integer.h>
#include <symengine/symengine_rcp.h>

// 引入枚举定义
#include "Region.hpp"
#include "Boundarys.hpp"

class RegionsInput
{
public:
    // 简单的矩形结构体
    struct Rect
    {
        double xl, xr, yb, yt;

        bool operator==(const Rect &other) const
        {
            return std::abs(xl - other.xl) < 1e-9 && std::abs(xr - other.xr) < 1e-9 &&
                   std::abs(yb - other.yb) < 1e-9 && std::abs(yt - other.yt) < 1e-9;
        }
    };

    struct SpecialRegionData
    {
        Rect area;
        double mu_r;
        double J_r;
    };

    // --- Properties ---
    std::vector<SpecialRegionData> special_regions; // 对应 special_region_area 和 property
    double air_mu_r;

    Rect all_area; // 所有需要计算的区域坐标

    std::vector<Rect> divided_rects; // 划分好的区域

    std::vector<CaseType> all_casetype;

    // 邻接区域索引 (0-based internally, output can be adjusted)
    std::vector<std::vector<int>> all_lefts;
    std::vector<std::vector<int>> all_rights;
    std::vector<std::vector<int>> all_tops;
    std::vector<std::vector<int>> all_bottoms;

    std::vector<int> all_BC_types; // 存储 BC_TYPE 的整数值或枚举

    int regions_num;
    int special_regions_num;

    int H_max = 60;
    int N_max = 60;

    std::vector<int> all_H_max;
    std::vector<int> all_N_max;
    std::vector<double> all_mu_r_val; // 存储数值方便计算
    std::vector<double> all_J_r_val;

    // --- Methods ---
    RegionsInput();
    virtual ~RegionsInput() = default;

    void set_current_regions(double xl, double xr, double yb, double yt, double mu_r, double J_r);

    void set_calculate_area(double xl, double xr, double yb, double yt, double mu_r);

    void divide_regions();

    void findNeighbors();

    void cal_other_info();

    // --- Return Functions (Getters) ---

    // 返回 SymEngine 格式的 regions (rects)
    std::vector<std::vector<SymEngine::RCP<const SymEngine::Symbol>>> rtn_regions_area();

    int rtn_region_num();

    // 返回特殊区域列表 (double 格式)
    std::vector<std::vector<double>> rtn_current_regions();

    std::pair<std::vector<int>, std::vector<int>> rtn_HN_max();

    // 返回 SymEngine 类型的 mu_r 和 J_r 以适配 AllRegions
    std::pair<std::vector<SymEngine::RCP<const SymEngine::Symbol>>,
              std::vector<SymEngine::RCP<const SymEngine::Basic>>>
    rtn_mu_J();

    std::pair<std::vector<int>, std::vector<CaseType>> rtn_types();

    // 返回边界索引，start_idx 用于偏移 (例如处理 yoz 平面时的索引偏移)
    struct BoundaryIndices
    {
        std::vector<std::vector<int>> lefts;
        std::vector<std::vector<int>> rights;
        std::vector<std::vector<int>> tops;
        std::vector<std::vector<int>> bottoms;
    };
    BoundaryIndices rtn_boundarys_idx(int start_idx);

    std::vector<int> rtn_current_idx(int start_idx);

private:
    bool intervalsOverlap(const std::pair<double, double> &a, const std::pair<double, double> &b);

    // 辅助: 浮点数比较
    bool eq(double a, double b) { return std::abs(a - b) < 1e-12; }
};