#pragma once

#include <vector>
#include <map>
#include <memory>
#include <string>

// SymEngine
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

#include "EnumTypes.h"

// Forward Declaration
class Region;
class BasicCase;
class Case1;
class Case2;

class Boundarys
{
public:
    // --- Properties ---
    std::vector<int> top;
    std::vector<int> bottom;
    std::vector<int> left;
    std::vector<int> right;

    BC_TYPE bc_type;

    // 指向拥有该 Boundarys 的 Region (避免循环引用，建议使用 weak_ptr 或原始指针)
    // MATLAB: obj.impl = impl (Region object)
    Region *impl;

    // 指向 Region 下的具体 Case 实现 (BasicCase)
    // MATLAB: obj.case_impl = impl.impl
    std::shared_ptr<BasicCase> case_impl;

    // 方程位置记录
    // MATLAB: BCfuncs_loc = zeros(2, all_region_num)
    // Key: region_idx, Value: pair(row, col) or just vector of pairs
    std::map<int, std::pair<int, int>> BCfuncs_loc;

    std::map<int, bool> BCfuncs_loc_map; // idx -> is_cal ?

    int all_region_num;

    // 区域到边界的映射
    // MATLAB: region2edge (索引->值)
    // Using map for sparse storage
    std::map<int, int> region2edge;

    // --- Methods ---
    Boundarys(Region *impl, int all_region_num);
    virtual ~Boundarys() = default;

    void cal_BC(bool Ln, bool Rn, bool Tn, bool Bn);
    void cal_ES(bool Ln, bool Rn, bool Tn, bool Bn);

private:
    // 辅助函数: genAA
    // 返回值: pair<coeffs list, equations list>
    struct GenResult
    {
        std::vector<SymEngine::RCP<const SymEngine::Basic>> coeffs;
        std::vector<SymEngine::RCP<const SymEngine::Basic>> eq;
    };

    GenResult genAA(int right_idx);

    // 辅助函数: genBB
    GenResult genBB(int right_idx, int cd_or_ef);

    // 辅助函数: getBB_func
    SymEngine::RCP<const SymEngine::Basic> getBB_func(const std::shared_ptr<BasicCase> &edge_impl, int x_or_y, int cd_or_ef);
};