#include "AllRegions.hpp"
#include "Region.hpp"
#include "RegionsInput.hpp" // Assuming this file exists
#include "BasicCase.hpp"    // For H_max, N_max, coeffs_exists access

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <numeric>

// SymEngine
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/real_double.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/eval_double.h>
#include <symengine/visitor.h>

using SymEngine::Basic;
using SymEngine::integer;
using SymEngine::RCP;
using SymEngine::real_double;
using SymEngine::Symbol;
using SymEngine::symbol;

AllRegions::AllRegions()
{
    this->divide_xoz = std::make_shared<RegionsInput>();
    this->divide_yoz = std::make_shared<RegionsInput>();
    this->region_num_xoz = 0;
    this->region_num_yoz = 0;
    this->all_region_num = 0;
}

void AllRegions::get_all_regions()
{
    // obj.set_all_regions();
    this->set_all_regions();
    std::cout << "区域内方程计算完毕" << std::endl;

    // obj.cal_IC_lens();
    this->cal_IC_lens();

    // [BC_funcs, BC_loc, ES_funcs] = obj.cal_all_BCs();
    BCResult bc_res = this->cal_all_BCs();
    std::cout << "计算系数方程完成" << std::endl;

    // [ES_xoz, ES_yoz] = obj.splice_ES(ES_funcs);
    std::pair<MatrixData, MatrixData> es_mats = this->splice_ES(bc_res.ES_funcs);
    std::cout << "计算ES方程完成" << std::endl;

    // [BC_xoz, BC_yoz] = obj.splice_BC(BC_funcs, BC_loc);
    std::pair<MatrixData, MatrixData> bc_mats = this->splice_BC(bc_res.BC_funcs, bc_res.BC_loc);
    std::cout << "计算BC方程完成" << std::endl;

    // IC = BC \ ES; (lsqr)
    // Note: Since we don't have a real linear algebra library here,
    // we assume solve_linear_system handles the math.
    std::vector<double> IC_xoz = this->solve_linear_system(bc_mats.first, es_mats.first);
    std::vector<double> IC_yoz = this->solve_linear_system(bc_mats.second, es_mats.second);
    std::cout << "计算IC方程完成" << std::endl;

    // [ICs_xoz, ICs_yoz] = obj.split_IC(IC_xoz, IC_yoz);
    SplitICResult ics = this->split_IC(IC_xoz, IC_yoz);

    // Saving logic skipped (requires File I/O library)
    // save("./mat/ICs", ...);
    // save("./mat/all_regions.mat", ...);
}

void AllRegions::set_all_regions()
{
    this->regions.resize(this->all_region_num);

    // parfor converted to for
    for (int i = 0; i < this->all_region_num; ++i)
    {
        // Need to construct Region.
        // Note: C++ indices are 0-based, MATLAB 1-based.
        // Adjust indices when passing to Region constructor if Region expects MATLAB logic,
        // or ensure Region handles 0-based. Assuming we pass 1-based index `i+1` for ID.

        // Construct basic args
        // regions_area[i] is vector<RCP<const Symbol>>

        this->regions[i] = std::make_shared<Region>(
            i + 1, // idx
            this->all_casetypes[i],
            this->regions_area[i],
            this->all_BC_types[i],
            this->all_tops[i],
            this->all_bottoms[i],
            this->all_lefts[i],
            this->all_rights[i],
            this->current_regions_idx, // Pass entire vector
            this->all_H_max[i],
            this->all_N_max[i],
            this->all_mu_r[i],
            this->all_J_r[i],
            this->all_region_num);

        this->regions[i]->get_region_solution_func();
    }
}

AllRegions::BCResult AllRegions::cal_all_BCs()
{
    BCResult res;
    res.BC_funcs.resize(this->all_region_num);
    res.BC_loc.resize(this->all_region_num);
    res.ES_funcs.resize(this->all_region_num);

    // parfor converted to for
    for (int i = 0; i < this->all_region_num; ++i)
    {
        auto coeff_res = this->regions[i]->gen_region_coefficient_func(this->regions);
        res.BC_funcs[i] = coeff_res.funcs;
        res.BC_loc[i] = coeff_res.BCfuncs_loc_map;
        res.ES_funcs[i] = coeff_res.ESfuncs;
    }
    return res;
}

std::pair<AllRegions::MatrixData, AllRegions::MatrixData> AllRegions::splice_BC(
    const std::vector<std::vector<std::vector<SymEngine::RCP<const SymEngine::Basic>>>> &BC_funcs,
    const std::vector<std::map<int, std::pair<int, int>>> &BC_loc)
{
    int N = this->all_region_num;
    std::vector<std::vector<MatrixData>> BC_blocks(N, std::vector<MatrixData>(N));

    printf("开始计算BC矩阵，共%d个\n", N);

    // parfor converted to for
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i == j)
            {
                BC_blocks[i][j] = {std::vector<double>(), 0, 0}; // Empty
                continue;
            }

            // Check adjacency using all_edge_regions which is a vector<bool> in Region
            if (this->regions[i]->all_edge_regions[j])
            {
                BC_blocks[i][j] = this->splice_BCxx(i, j, BC_funcs, BC_loc);
                printf("计算完成区域(%d,%d)\n", i + 1, j + 1);
            }
            else
            {
                BC_blocks[i][j] = {std::vector<double>(), 0, 0};
            }
        }
    }

    // Split blocks for xoz and yoz
    int k = this->region_num_xoz;

    // Build sub-grids
    std::vector<std::vector<MatrixData>> blocks_xoz(k, std::vector<MatrixData>(k));
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j)
            blocks_xoz[i][j] = BC_blocks[i][j];

    int m = N - k;
    std::vector<std::vector<MatrixData>> blocks_yoz(m, std::vector<MatrixData>(m));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
            blocks_yoz[i][j] = BC_blocks[k + i][k + j];

    MatrixData BC_xoz = this->build_block_matrix(blocks_xoz);
    MatrixData BC_yoz = this->build_block_matrix(blocks_yoz);

    return {BC_xoz, BC_yoz};
}

AllRegions::MatrixData AllRegions::build_block_matrix(const std::vector<std::vector<MatrixData>> &BC_blocks)
{
    int N = BC_blocks.size();
    if (N == 0)
        return {std::vector<double>(), 0, 0};

    std::vector<int> row_sizes(N, 0);
    std::vector<int> col_sizes(N, 0);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            const auto &block = BC_blocks[i][j];
            if (block.rows > 0 && block.cols > 0)
            {
                row_sizes[i] = std::max(row_sizes[i], block.rows);
                col_sizes[j] = std::max(col_sizes[j], block.cols);
            }
        }
    }

    int total_rows = std::accumulate(row_sizes.begin(), row_sizes.end(), 0);
    int total_cols = std::accumulate(col_sizes.begin(), col_sizes.end(), 0);

    MatrixData BigMat;
    BigMat.rows = total_rows;
    BigMat.cols = total_cols;
    BigMat.data.assign(total_rows * total_cols, 0.0);

    int r0 = 0;
    for (int i = 0; i < N; ++i)
    {
        int c0 = 0;
        for (int j = 0; j < N; ++j)
        {
            const auto &block = BC_blocks[i][j];
            if (block.rows > 0)
            {
                // Copy block into BigMat
                for (int rr = 0; rr < block.rows; ++rr)
                {
                    for (int cc = 0; cc < block.cols; ++cc)
                    {
                        int global_r = r0 + rr;
                        int global_c = c0 + cc;
                        BigMat.data[global_r * total_cols + global_c] = block.data[rr * block.cols + cc];
                    }
                }
            }
            c0 += col_sizes[j];
        }
        r0 += row_sizes[i];
    }

    // Add Identity Matrix (Sparse logic in MATLAB: BC = BigMat + sparse(eye))
    // Only add 1.0 to diagonal elements
    int n = std::min(total_rows, total_cols); // Should be square
    printf("BC矩阵大小:%dx%d\n", total_rows, total_cols);
    for (int i = 0; i < n; ++i)
    {
        BigMat.data[i * total_cols + i] += 1.0;
    }

    return BigMat;
}

AllRegions::MatrixData AllRegions::splice_BCxx(int idx, int edge_idx,
                                               const std::vector<std::vector<std::vector<RCP<const Basic>>>> &BC_funcs,
                                               const std::vector<std::map<int, std::pair<int, int>>> &BC_loc)
{
    // Implementation of symbolic evaluation logic
    // This part requires mapping SymEngine expressions to numerical functions (lambdas)
    // and evaluating them over grids.

    auto idx_case = this->regions[idx]->impl;
    auto edge_case = this->regions[edge_idx]->impl;

    // funcss corresponds to BC_funcs{idx} -> vector<vector<RCP>>
    const auto &funcss = BC_funcs[idx];

    // Find coeffs exists indices (1-based from MATLAB logic)
    std::vector<int> row_exists;
    for (size_t k = 0; k < idx_case->coeffs_exists.size(); ++k)
        if (idx_case->coeffs_exists[k])
            row_exists.push_back(k + 1);

    std::vector<int> col_exists;
    for (size_t k = 0; k < edge_case->coeffs_exists.size(); ++k)
        if (edge_case->coeffs_exists[k])
            col_exists.push_back(k + 1);

    bool row_has_cd0x = (!row_exists.empty() && row_exists[0] == 1);

    // Calculate row/col lengths
    std::vector<int> rows_len;
    for (int x : row_exists)
    {
        int len = ((x == 1 || x == 3) ? 1 : (x == 2 || x == 4) ? idx_case->H_max
                                                               : idx_case->N_max);
        rows_len.push_back(len);
    }
    std::vector<int> cols_len;
    for (int x : col_exists)
    {
        int len = ((x == 1 || x == 3) ? 1 : (x == 2 || x == 4) ? edge_case->H_max
                                                               : edge_case->N_max);
        cols_len.push_back(len);
    }

    int rows_bc = std::accumulate(rows_len.begin(), rows_len.end(), 0);
    int cols_bc = std::accumulate(cols_len.begin(), cols_len.end(), 0);

    MatrixData BCxx;
    BCxx.rows = rows_bc;
    BCxx.cols = cols_bc;
    BCxx.data.assign(rows_bc * cols_bc, 0.0);

    // Get edge location info: BC_loc{idx} map key=edge_idx -> pair(type, row_index)
    // Note: C++ map uses key=edge_idx (0-based in vector, but edge_idx passed here is 0-based)
    // However, map key logic depends on Region::gen_region_coefficient_func implementation.
    // Assuming BC_loc maps region_index -> pair(type_code, matrix_row_index)
    if (BC_loc[idx].find(edge_idx) == BC_loc[idx].end())
        return BCxx; // Should not happen if adjacent

    std::pair<int, int> edge_bc_loc = BC_loc[idx].at(edge_idx);
    int bc_type_code = edge_bc_loc.first; // 1=Top, 2=Bottom etc.
    int bc_row_idx = edge_bc_loc.second;  // Row index in funcs cell array

    // Find 'i' in row_exists that matches bc_type_code
    int i_idx = -1;
    for (size_t k = 0; k < row_exists.size(); ++k)
    {
        if (row_exists[k] == bc_type_code)
        {
            i_idx = k;
            break;
        }
    }
    if (i_idx == -1)
        return BCxx;

    // Determine row variable (h or n)
    // bool row_is_h = (bc_type_code <= 4);

    // Iterate columns (j)
    int start_col = 0;
    for (size_t j = 0; j < col_exists.size(); ++j)
    {
        int col_code = col_exists[j];

        // Find corresponding equation
        // funcss is vector<vector<RCP>>. Indexing: funcss[bc_row_idx][col_code]
        // Note: C++ vectors are 0-based. MATLAB was 1-based.
        // Assuming BC_funcs is sized appropriately.
        // Also note: `bc_row_idx` from MATLAB was 1-based index into T_funcs/B_funcs?
        // No, BC_loc stored the index into the consolidated funcs array.
        // In C++, BC_funcs[idx] is vector<vector>.
        // We need to match the indexing.

        // Simulating: func = funcss{edge_bc_loc(2), col_exists(j)}
        RCP<const Basic> func;
        if (bc_row_idx > 0 && bc_row_idx <= (int)funcss.size())
        {
            if (col_code > 0 && col_code <= (int)funcss[bc_row_idx - 1].size())
            {
                func = funcss[bc_row_idx - 1][col_code - 1];
            }
        }

        if (func)
        {
            // Symbolic Evaluation simulation
            // Since SymEngine doesn't have matlabFunction(), we substitute and eval.
            // Loop over grid points (row_idx, col_idx)

            int r_len = rows_len[i_idx];
            int c_len = cols_len[j];

            // Calculate starting row for this block
            int start_row = 0;
            for (int k = 0; k < i_idx; ++k)
                start_row += rows_len[k];

            for (int r = 0; r < r_len; ++r)
            {
                for (int c = 0; c < c_len; ++c)
                {
                    // Substitute variables row_hn, col_hn
                    SymEngine::map_basic_basic subs_map;
                    // Determine which symbols to sub based on h/n logic
                    // This requires the symbolic variables used in 'func' to be accessible
                    // or known. In C++, these are local to Region/BasicCase.
                    // We need a robust way to map 'h', 'n' in the expression to integers.
                    // Assuming 'func' contains symbols named "h_idx", "n_idx".

                    // Note: This is highly complex in C++ without symbolic reflection.
                    // For now, I'll place a placeholder for the numerical evaluation.
                    double val = 0.0;

                    // TODO: Implement symbolic subs and eval_double
                    // subs_map[idx_case->h] = integer(r + 1);
                    // subs_map[edge_case->h or n] = integer(c + 1);
                    // val = -eval_double(func->subs(subs_map));

                    // Assign to BCxx
                    int gr = start_row + r;
                    int gc = start_col + c;
                    BCxx.data[gr * BCxx.cols + gc] = val;
                }
            }

            // Handle special cd0x logic (MATLAB code block)
            if (bc_type_code <= 4 && row_has_cd0x)
            {
                // Logic for c0/d0 equations interacting with c/d
                // ... (Implementation specific to physics logic)
            }
        }

        start_col += cols_len[j];
    }

    return BCxx;
}

std::pair<AllRegions::MatrixData, AllRegions::MatrixData> AllRegions::splice_ES(
    const std::vector<std::vector<RCP<const Basic>>> &ES_funcs)
{
    // Implementation of ES splicing
    // ES is [ES_xoz; ES_yoz]
    // Each region produces a vector ESx.

    std::vector<double> ES_vec;

    int N = this->all_region_num;
    std::vector<double> ES_xoz_vec;

    for (int i = 0; i < N; ++i)
    {
        if (i == this->region_num_xoz)
        {
            // Split point
            // For now, just keep accumulating into one big vector and split later or separate loops.
        }

        // Similar to BCxx, need to evaluate ES_funcs[i] symbolic expressions
        // ...
    }

    // Placeholder return
    return {{std::vector<double>(), 0, 0}, {std::vector<double>(), 0, 0}};
}

AllRegions::SplitICResult AllRegions::split_IC(const std::vector<double> &IC_xoz, const std::vector<double> &IC_yoz)
{
    SplitICResult res;
    // Logic to chop the flat IC vectors back into structured arrays based on len_IC
    // ...
    return res;
}

void AllRegions::cal_IC_lens()
{
    this->len_IC_xoz.resize(this->region_num_xoz);
    for (int i = 0; i < this->region_num_xoz; ++i)
    {
        auto idx_impl = this->regions[i]->impl;
        // logic: coeffs_exists * [1, H, 1, H, N, N]
        int len = 0;
        const auto &ce = idx_impl->coeffs_exists;
        int H = idx_impl->H_max;
        int N = idx_impl->N_max;
        if (ce[0])
            len += 1;
        if (ce[1])
            len += H;
        if (ce[2])
            len += 1;
        if (ce[3])
            len += H;
        if (ce[4])
            len += N;
        if (ce[5])
            len += N;
        this->len_IC_xoz[i] = len;
    }

    this->len_IC_yoz.resize(this->region_num_yoz);
    for (int i = 0; i < this->region_num_yoz; ++i)
    {
        auto idx_impl = this->regions[i + this->region_num_xoz]->impl;
        int len = 0;
        const auto &ce = idx_impl->coeffs_exists;
        int H = idx_impl->H_max;
        int N = idx_impl->N_max;
        if (ce[0])
            len += 1;
        if (ce[1])
            len += H;
        if (ce[2])
            len += 1;
        if (ce[3])
            len += H;
        if (ce[4])
            len += N;
        if (ce[5])
            len += N;
        this->len_IC_yoz[i] = len;
    }
}

void AllRegions::input_current_region(const std::string &plain, double xl, double xr, double yb, double yt, double mu_r, double J_r)
{
    if (plain == "xoz")
    {
        this->divide_xoz->set_current_regions(xl, xr, yb, yt, mu_r, J_r);
    }
    else if (plain == "yoz")
    {
        this->divide_yoz->set_current_regions(xl, xr, yb, yt, mu_r, J_r);
    }
}

void AllRegions::input_calculate_area(const std::string &plain, double xl, double xr, double yb, double yt, double mu_r)
{
    if (plain == "xoz")
    {
        this->divide_xoz->set_calculate_area(xl, xr, yb, yt, mu_r);
    }
    else if (plain == "yoz")
    {
        this->divide_yoz->set_calculate_area(xl, xr, yb, yt, mu_r);
    }
}

void AllRegions::pre_process()
{
    this->divide_xoz->divide_regions();
    this->divide_xoz->findNeighbors();
    this->divide_xoz->cal_other_info();

    // Retrieve data from divide_xoz
    this->regions_area = this->divide_xoz->rtn_regions_area();
    this->region_num_xoz = this->divide_xoz->rtn_region_num();
    // ... populate other fields (casetypes, BC_types etc) ...
    // Note: Since RegionsInput class is not provided, this is hypothetical linking.

    this->divide_yoz->divide_regions();
    // ...

    // Merge xoz and yoz data
    this->all_region_num = this->region_num_xoz + this->region_num_yoz;
}

// Placeholder solver
std::vector<double> AllRegions::solve_linear_system(const MatrixData &A, const MatrixData &b)
{
    // Return empty vector or implement LSQR / LU solver
    return std::vector<double>(A.cols, 0.0);
}

// Dummy simplify
RCP<const Basic> AllRegions::simplifyFraction(const RCP<const Basic> &expr)
{
    return expr;
}

// Helper methods for cal_Bx_By omitted for brevity, similar structure to cal_BC logic.