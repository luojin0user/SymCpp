#pragma once

#include <vector>
#include <memory>
#include <string>
#include <stdexcept>

// 引入 SymEngine
#include <symengine/symbol.h>
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>

// 引入所有具体的 Case 实现
#include "BasicCase.hpp"
#include "Case1.hpp"
#include "Case2.hpp"

#include "EnumTypes.h"
#include "Boundarys.hpp"

#include "RegionsInput.hpp"

class Region
{
private:
    int boundary_status; // 存储边界类型

    // 边界管理对象
    std::unique_ptr<Boundarys> boundarys;

    int all_regions_num;

    void set_boundary(const std::vector<int> &top,
                      const std::vector<int> &bottom,
                      const std::vector<int> &left,
                      const std::vector<int> &right,
                      BC_TYPE bc_type);

public:
    CaseType case_type;

    // 储存所有区域的节点
    const std::vector<std::unique_ptr<Region>> *all_regions;

    // 所有邻接区域的数组
    std::vector<std::pair<bool, int>> all_edge_regions;

    // 多态指针，指向具体实现 (BasicCase 子类)
    std::unique_ptr<BasicCase> impl;

    // Q矩阵最后处在BCxx矩阵中的位置（大的行），记录的当前idx的所有邻接的edge(c0 c d0 d e f)在对应的矩阵的位置（值）
    int BCxx_loc[6] = {-1, -1, -1, -1, -1, -1};

    // --- Methods ---
    Region(int idx,
           CaseType type,
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
           int all_regions_num);

    virtual ~Region() = default;

    // 返回方程组 {eq_A_z, eq_B_x}
    void get_region_solution_func();

    // 生成系数方程
    void gen_region_coefficient_func(const std::vector<std::unique_ptr<Region>> *all_regions);

    void cal_BCxx_loc();
};