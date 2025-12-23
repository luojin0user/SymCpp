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
#include "EnumTypes.h"

class RegionsInput
{
public:
    struct SpecialRegionData
    {
        struct Rect area;
        double mu_r;
        double J_r;
    };

    // --- Properties ---
    std::vector<Rect> special_regions_area; // 对应 special_region_area 和 property
    std::vector<std::pair<double, double>> special_regions_property;

    struct Rect all_area; // 所有需要计算的区域坐标
    double air_mu_r;

    std::vector<Rect> divided_rects; // 划分好的区域

    std::vector<CaseType> all_casetype;

    // 邻接区域索引 (0-based internally, output can be adjusted)
    std::vector<std::vector<int>> all_lefts;
    std::vector<std::vector<int>> all_rights;
    std::vector<std::vector<int>> all_tops;
    std::vector<std::vector<int>> all_bottoms;

    std::vector<BC_TYPE> all_BC_types; // 存储 BC_TYPE 枚举

    int regions_num;
    int special_regions_num;

    int H_max = 1;
    int N_max = 1;

    std::vector<int> all_H_max;
    std::vector<int> all_N_max;
    std::vector<float> all_mu_r_val; // 存储数值方便计算
    std::vector<float> all_J_r_val;

    std::vector<int> current_idx;

    // --- Methods ---
    RegionsInput();
    virtual ~RegionsInput() = default;

    void set_current_regions(float xl, float xr, float yb, float yt, float mu_r, float J_r);

    void set_calculate_area(float xl, float xr, float yb, float yt, float mu_r);

    void divide_regions();

    void findNeighbors();

    void cal_other_info();

    // --- Return Functions (Getters) ---
    void set_boundarys_idx(int start_idx);
    void cal_current_idx(int start_idx);

    // 返回该平面所有的信息
    struct Plane_Info rtn_Plane_Info(int start_idx);

private:
    bool intervalsOverlap(const std::pair<float, float> &a, const std::pair<float, float> &b);

    // 辅助: 浮点数比较
    bool eq(float a, float b) { return std::abs(a - b) < 1e-8; }
};