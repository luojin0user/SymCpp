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

class Region
{
public:
    // --- Properties ---
    int idx;
    SymEngine::RCP<const SymEngine::Symbol> mu_r;
    SymEngine::RCP<const SymEngine::Basic> J_r;

    CaseType case_type;
    int boundary_status; // 存储边界类型

    // 多态指针，指向具体实现 (BasicCase 子类)
    std::shared_ptr<BasicCase> impl;

    // 边界管理对象
    std::shared_ptr<Boundarys> boundarys;

    // 储存所有区域的节点 (使用 shared_ptr 避免拷贝)
    std::vector<std::shared_ptr<Region>> all_regions;

    // 边界存在性标志 (True 代表该边界不存在/是开放的?, 根据 MATLAB注释: "如果为True，说明这个边界不存在")
    // 实际上通常 Ln=True 意味着 IsNeighbor (有邻居)，所以该物理边界不存在。
    bool Ln, Rn, Tn, Bn;

    int all_regions_num;

    // 所有邻接区域的数组
    std::vector<bool> all_edge_regions;

    // --- Methods ---
    Region(int idx, CaseType type,
           const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &area,
           int bc_type,
           const std::vector<int> &top,
           const std::vector<int> &bottom,
           const std::vector<int> &left,
           const std::vector<int> &right,
           const std::vector<int> &ES_regions_idx,
           int H_max, int N_max,
           const SymEngine::RCP<const SymEngine::Symbol> &mu_r,
           const SymEngine::RCP<const SymEngine::Basic> &J_r,
           int all_regions_num);

    virtual ~Region() = default;

    void set_boundary(const std::vector<int> &top,
                      const std::vector<int> &bottom,
                      const std::vector<int> &left,
                      const std::vector<int> &right,
                      int bc_type);

    // 返回方程组 {eq_A_z, eq_B_x}
    std::vector<SymEngine::RCP<const SymEngine::Basic>> get_region_solution_func();

    // 生成系数方程
    // 返回值: tuple<系数方程列表, 位置映射, 源项方程>
    // funcs: 这里为了通用性，将 cell 数组展平或按组返回
    struct CoefficientResult
    {
        std::vector<std::vector<SymEngine::RCP<const SymEngine::Basic>>> funcs;
        std::map<int, std::pair<int, int>> BCfuncs_loc_map;
        std::vector<SymEngine::RCP<const SymEngine::Basic>> ESfuncs;
    };

    CoefficientResult gen_region_coefficient_func(const std::vector<std::shared_ptr<Region>> &all_regions);
};