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

    // 指向拥有该 Boundarys 的 Region
    const Region *region;
    const std::unique_ptr<BasicCase> &region_impl;

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
    Boundarys(const Region *impl, int all_region_num);
    virtual ~Boundarys() = default;

    void cal_BC(bool Ln, bool Rn, bool Tn, bool Bn);
    void cal_ES(bool Ln, bool Rn, bool Tn, bool Bn);

    // 将这个类中的数据写入到外界的容器中
    void rtn_BC_funcs(std::vector<Boundary_Funcs> &T_funcs,
                      std::vector<Boundary_Funcs> &B_funcs,
                      std::vector<Boundary_Funcs> &L_funcs,
                      std::vector<Boundary_Funcs> &R_funcs);

    void rtn_ES_funcs(std::vector<SymEngine::RCP<const SymEngine::Basic>> &T_ESfuncs,
                      std::vector<SymEngine::RCP<const SymEngine::Basic>> &B_ESfuncs,
                      std::vector<SymEngine::RCP<const SymEngine::Basic>> &L_ESfuncs,
                      std::vector<SymEngine::RCP<const SymEngine::Basic>> &R_ESfuncs);

private:
    // 辅助函数: genAA
    // 返回值: pair<coeffs list, equations list>
    struct GenResult
    {
        std::vector<SymEngine::RCP<const SymEngine::Basic>> coeffs;
        std::vector<SymEngine::RCP<const SymEngine::Basic>> eq;
    };

    void process_boundary_BB(
        bool BC_exist,                            // 边界是否存在 !Bn，如果为true则存在边界
        const std::vector<int> &boundary_indices, // 当前处理的的边界
        const SymEngine::RCP<const SymEngine::Basic> &sub_src,
        const SymEngine::RCP<const SymEngine::Basic> &sub_dst,
        std::vector<Boundary_Funcs> &funcs,                             // 需要返回的边界方程
        std::vector<SymEngine::RCP<const SymEngine::Basic>> &ES_funcs); // 可能存在的ES方程

    void process_boundary_AA(
        bool BC_exist,                            // 边界是否存在 !Bn，如果为true则存在边界
        const std::vector<int> &boundary_indices, // 当前处理的的边界
        const SymEngine::RCP<const SymEngine::Basic> &sub_src,
        const SymEngine::RCP<const SymEngine::Basic> &sub_dst,
        int tb_or_lr,                                                   // 当前处理的是上下(tb,1)或者左右(lr,2)
        std::vector<Boundary_Funcs> &funcs,                             // 需要返回的边界方程
        std::vector<SymEngine::RCP<const SymEngine::Basic>> &ES_funcs); // 可能存在的ES方程

    Boundary_Funcs genAA(int right_idx);

    // 辅助函数: genBB
    Boundary_Funcs genBB(int right_idx, int cd_or_ef);

    // 辅助函数: getBB_func
    SymEngine::RCP<const SymEngine::Basic> getBB_func(const BasicCase &edge_impl, int x_or_y, int cd_or_ef);

    // 边界函数
    std::vector<Boundary_Funcs> B_funcs;
    std::vector<Boundary_Funcs> T_funcs;
    std::vector<Boundary_Funcs> L_funcs;
    std::vector<Boundary_Funcs> R_funcs;

    // 有源项函数 (ES_funcs)
    std::vector<SymEngine::RCP<const SymEngine::Basic>> T_ESfuncs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> B_ESfuncs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> L_ESfuncs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> R_ESfuncs;
};