#pragma once

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <algorithm>

// SymEngine
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

// Forward Declarations
class Region;
class RegionsInput; // Assuming this class exists based on MATLAB code
enum class CaseType;
enum class BC_TYPE;

class AllRegions
{
public:
    // --- Properties ---

    // 存储所有区域 (Use shared_ptr for automatic memory management)
    std::vector<std::shared_ptr<Region>> regions;

    // 用于处理输入区域等 (Assuming RegionsInput is a defined class)
    std::shared_ptr<RegionsInput> divide_xoz;
    std::shared_ptr<RegionsInput> divide_yoz;

    int region_num_xoz; // 区域数量
    int region_num_yoz; // 区域数量
    int all_region_num;

    std::vector<int> all_H_max;
    std::vector<int> all_N_max;

    std::vector<SymEngine::RCP<const SymEngine::Symbol>> all_mu_r;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> all_J_r;

    std::vector<int> len_IC_xoz;
    std::vector<int> len_IC_yoz;

    // 区域，每一行是一个矩形区域 [xl, xr, yl, yt]
    // Mapped to vector of vector<Symbol> or double depending on usage.
    // Based on previous files, Region constructor takes vector<RCP<const Symbol>>
    std::vector<std::vector<SymEngine::RCP<const SymEngine::Symbol>>> regions_area;

    std::vector<int> all_BC_types; // Storing int representation of BC_TYPE
    std::vector<CaseType> all_casetypes;

    // 有源区域
    // MATLAB: current_regions (Unclear type, possibly bool array or list of indices)
    // Based on usage [current_regions1; current_regions2], likely a list of properties or boolean flags.
    // Let's assume boolean mask for now given usage in splice_ES?
    // Or actually, MATLAB code has `obj.current_regions_idx`, so `current_regions` might be the raw input data.
    // I will use vector<int> for simplicity or vector<bool> if it's a mask.
    std::vector<int> current_regions;
    std::vector<int> current_regions_idx;

    // 所有区域的上下左右边界 (Adjacency lists)
    std::vector<std::vector<int>> all_lefts;
    std::vector<std::vector<int>> all_rights;
    std::vector<std::vector<int>> all_tops;
    std::vector<std::vector<int>> all_bottoms;

    // 当前平面，xoy或者zoy字符串
    std::string this_plain;

    // --- Methods ---

    AllRegions();
    virtual ~AllRegions() = default;

    void get_all_regions();

    void set_all_regions();

    // Return tuple: BC_funcs, BC_loc, ES_funcs
    struct BCResult
    {
        std::vector<std::vector<std::vector<SymEngine::RCP<const SymEngine::Basic>>>> BC_funcs;
        std::vector<std::map<int, std::pair<int, int>>> BC_loc;
        std::vector<std::vector<SymEngine::RCP<const SymEngine::Basic>>> ES_funcs;
    };

    BCResult cal_all_BCs();

    // 拼接所有的BC
    // Returns pair of sparse matrices (represented as dense for now or custom struct)
    // Since we don't have a sparse matrix library here, we return a flat representation or void
    // Here I will use a placeholder double matrix structure for the solver.
    struct MatrixData
    {
        // Placeholder for a dense matrix (row-major)
        std::vector<double> data;
        int rows;
        int cols;
    };

    std::pair<MatrixData, MatrixData> splice_BC(
        const std::vector<std::vector<std::vector<SymEngine::RCP<const SymEngine::Basic>>>> &BC_funcs,
        const std::vector<std::map<int, std::pair<int, int>>> &BC_loc);

    MatrixData build_block_matrix(const std::vector<std::vector<MatrixData>> &BC_blocks);

    MatrixData splice_BCxx(int idx, int edge_idx,
                           const std::vector<std::vector<std::vector<SymEngine::RCP<const SymEngine::Basic>>>> &BC_funcs,
                           const std::vector<std::map<int, std::pair<int, int>>> &BC_loc);

    // 拼接ES矩阵
    std::pair<MatrixData, MatrixData> splice_ES(const std::vector<std::vector<SymEngine::RCP<const SymEngine::Basic>>> &ES_funcs);

    // 分割IC矩阵
    struct SplitICResult
    {
        // Cell array equivalent: vector of vector of vectors (Region -> Coefficient Index -> Values)
        std::vector<std::vector<std::vector<double>>> ICs_xoz;
        std::vector<std::vector<std::vector<double>>> ICs_yoz;
    };

    SplitICResult split_IC(const std::vector<double> &IC_xoz, const std::vector<double> &IC_yoz);

    void cal_IC_lens();

    // Calculation Helpers
    void cal_Bx_By_region(const std::string &plane, const SplitICResult &ICs, int region_num,
                          const std::vector<double> &x0, const std::vector<double> &y0,
                          std::vector<double> &Bx, std::vector<double> &By);

    void cal_Bx_By(const std::string &plane, const SplitICResult &ICs,
                   const std::vector<double> &x0, const std::vector<double> &y0,
                   std::vector<double> &Bx_out, std::vector<double> &By_out);

    // Input proxies
    void input_current_region(const std::string &plain, double xl, double xr, double yb, double yt, double mu_r, double J_r);
    void input_calculate_area(const std::string &plain, double xl, double xr, double yb, double yt, double mu_r);

    void pre_process();

private:
    // Helper to solve Ax = b (Placeholder for lsqr)
    std::vector<double> solve_linear_system(const MatrixData &A, const MatrixData &b);

    // Helper for simplifyFraction (dummy implementation)
    SymEngine::RCP<const SymEngine::Basic> simplifyFraction(const SymEngine::RCP<const SymEngine::Basic> &expr);
};