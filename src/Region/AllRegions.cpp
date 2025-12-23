#include "AllRegions.hpp"
#include "Region.hpp"
#include "RegionsInput.hpp" // Assuming this file exists
#include "BasicCase.hpp"    // For H_max, N_max, coeffs_exists access
#include "Matrix2D.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <fstream>
#include <string>
#include <iomanip>

// SymEngine
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/real_double.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/eval_double.h>
#include <symengine/visitor.h>

#include <symengine/lambda_double.h>

using SymEngine::Basic;
using SymEngine::integer;
using SymEngine::RCP;
using SymEngine::real_double;
using SymEngine::Symbol;
using SymEngine::symbol;

AllRegions::AllRegions()
{
    this->divide_xoz = std::make_unique<RegionsInput>();
    this->divide_yoz = std::make_unique<RegionsInput>();
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

    int len_xoz = 0, len_yoz = 0;
    for (int i = 0; i < region_num_xoz; i++)
        len_xoz += all_len_IC[i];

    for (int i = region_num_xoz; i < region_num_xoz + region_num_yoz; i++)
        len_yoz += all_len_IC[i];

    // [BC_funcs, BC_loc, ES_funcs] = obj.cal_all_BCs();
    this->cal_all_BCs();
    std::cout << "计算系数方程完成" << std::endl;

    // [ES_xoz, ES_yoz] = obj.splice_ES(ES_funcs);
    Matrix2D ES_xoz(len_xoz, 1); // 已经全部初始化为0
    Matrix2D ES_yoz(len_yoz, 1);
    this->splice_ES(ES_xoz, ES_yoz);
    std::cout << "计算ES方程完成" << std::endl;

    // [BC_xoz, BC_yoz] = obj.splice_BC(BC_funcs, BC_loc);

    Matrix2D BC_xoz(len_xoz, len_xoz); // 已经全部初始化为0
    Matrix2D BC_yoz(len_yoz, len_yoz);
    this->splice_BC(BC_xoz, BC_yoz);
    std::cout << "计算BC方程完成" << std::endl;

    // IC = BC \ ES; (lsqr)
    auto ES1 = ES_xoz.toEigenCopy();
    auto ES2 = ES_yoz.toEigenCopy();
    auto BC1 = BC_xoz.toEigenCopy();
    auto BC2 = BC_yoz.toEigenCopy();

    auto IC1 = BC1.ldlt().solve(ES1);
    auto IC2 = BC2.ldlt().solve(ES2);

    saveEigenToCSV("BC1.csv", BC1);

    std::cout << "计算IC方程完成" << std::endl;

    // [ICs_xoz, ICs_yoz] = obj.split_IC(IC_xoz, IC_yoz);
    // SplitICResult ics = this->split_IC(IC_xoz, IC_yoz);

    // Saving logic skipped (requires File I/O library)
    // save("./mat/ICs", ...);
    // save("./mat/all_regions.mat", ...);
}

void AllRegions::set_all_regions()
{
    this->regions.resize(this->all_region_num);

    for (int i = 0; i < this->all_region_num; ++i)
    {
        this->regions[i] = std::make_unique<Region>(
            i, // idx
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

void AllRegions::cal_all_BCs()
{
    for (int i = 0; i < this->all_region_num; ++i)
    {
        this->regions[i]->gen_region_coefficient_func(&this->regions);
    }
}

void AllRegions::splice_BC(Matrix2D &BC_xoz, Matrix2D &BC_yoz)
{
    int N = this->all_region_num;

    std::vector<std::vector<Matrix2D>> BC_blocks_xoz;
    std::vector<std::vector<Matrix2D>> BC_blocks_yoz;
    BC_blocks_xoz.resize(region_num_xoz);
    BC_blocks_yoz.resize(region_num_yoz);

    for (auto &row : BC_blocks_xoz)
    {
        row.resize(region_num_xoz); // 默认构造，不拷贝
    }

    for (auto &row : BC_blocks_yoz)
    {
        row.resize(region_num_yoz); // 默认构造，不拷贝
    }

    printf("开始计算BC矩阵，共%d个\n", N);

#pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if ((i < region_num_xoz && j >= region_num_xoz) ||
                (i >= region_num_xoz && j < region_num_xoz))
                continue;
            if (i == j)
            {
                // BC_blocks[i][j] = {std::vector<std::vector<double>>(), 0, 0};
                continue;
            }

            // Check adjacency using all_edge_regions which is a vector<bool> in Region
            if (this->regions[i]->all_edge_regions[j].first)
            {
                // 首先计算这个数组有多少行，生成BCxx数组
                if (i < region_num_xoz)
                {
                    BC_blocks_xoz[i][j] = Matrix2D(this->all_len_IC[i], this->all_len_IC[j]);
                    this->splice_BCxx(i, j, BC_blocks_xoz[i][j]);
                }
                else
                {
                    BC_blocks_yoz[i - region_num_xoz][j - region_num_xoz] = Matrix2D(this->all_len_IC[i], this->all_len_IC[j]);
                    this->splice_BCxx(i, j, BC_blocks_yoz[i - region_num_xoz][j - region_num_xoz]);
                }

                printf("计算完成区域(%d,%d)\n", i, j);
            }
        }
    }

    // copy数据
    int row_s = 0, col_s = 0; // 当前开始的坐标

    for (int i = 0; i < region_num_xoz; ++i)
    {
        col_s = 0;
        for (int j = 0; j < region_num_xoz; ++j)
        {
            if (this->regions[i]->all_edge_regions[j].first)
            {
                BC_xoz.write_block(row_s, col_s, BC_blocks_xoz[i][j]);
            }
            col_s += all_len_IC[j];
        }

        row_s += all_len_IC[i];
    }

    row_s = 0, col_s = 0;
    for (int i = 0; i < region_num_yoz; ++i)
    {
        col_s = 0;
        for (int j = 0; j < region_num_yoz; ++j)
        {
            if (this->regions[region_num_xoz + i]->all_edge_regions[region_num_xoz + j].first)
            {
                BC_yoz.write_block(row_s, col_s, BC_blocks_yoz[i][j]);
            }

            col_s += all_len_IC[region_num_xoz + j];
        }

        row_s += all_len_IC[region_num_xoz + i];
    }
}

// 计算出一个BCxx
void AllRegions::splice_BCxx(int idx, int edge_idx, Matrix2D &BCxx)
{
    // 首先，直接调用region类的方程即可
    // std::array<std::vector<std::array<Integral_Func, 6>>, 6> eq_BC;
    // 最外层为上下左右边界 c0 c d0 d e f
    // 中间是这个边界的分段函数情况
    // 最内层是对应的邻接区域的c0 c d0 d e f

    // 为了知道2个邻接区域的方程，固然需要知道邻接区域与这个方程数组的映射关系
    // 也就是知道idx和edge_idx，映射到eq_BC的最外层的几行几列，最内层无需计算（这是单个区域的）
    // 对于每个idx都有自己的储存，唯一的变量就是edge_idx
    // std::vector<std::tuple<int, int, bool>> eq_BC_loc;
    const auto &idx_impl = regions[idx]->impl;
    const auto &edge_impl = regions[edge_idx]->impl;

    const auto &loc = idx_impl->eq_BC_loc[edge_idx];
    const auto &funcs = idx_impl->eq_BC[std::get<0>(loc)][std::get<1>(loc)];

    // 需要找到这个邻接区域在BCxx中的位置，即第几行开始

    // 然后确定h n的值，这个可以由loc确定row的，由funcs确定col的
    RCP<const Basic> row_hn;
    RCP<const Basic> col_hn;
    unsigned row_num;
    unsigned col_num; // 行数和列数

    if (std::get<0>(loc) <= 3)
    {
        row_hn = idx_impl->h;
        row_num = idx_impl->const_vals.H_max;
    }
    else
    {
        row_hn = idx_impl->n;
        row_num = idx_impl->const_vals.N_max;
    }

    // 循环计算对应的值
    for (int j = 0; j < funcs.size(); j++)
    {
        // 如果当前这个不存在，那么继续找下一个
        if (!edge_impl->coeffs_exists[j])
            continue;

        // 对面区域的h n及其数量
        if (j <= 3)
        {
            col_hn = edge_impl->h;
            col_num = edge_impl->const_vals.H_max;
        }
        else
        {
            col_hn = edge_impl->n;
            col_num = edge_impl->const_vals.N_max;
        }

        // 循环计算
        auto expr = funcs[j].funcs;
        double a = funcs[j].start_int;
        double b = funcs[j].end_int;

        // 转成std::function
        SymEngine::vec_basic inputs;
        inputs.push_back(funcs[j].x);

        // 如果这个是c0或者d0
        if (j == 0 && j == 2)
        {
            inputs.push_back(row_hn);
            col_num = 1;
        }
        else
        {
            inputs.push_back(row_hn);
            inputs.push_back(col_hn);
        }

        SymEngine::LambdaDoubleVisitor<long double> visitor;
        // std::cout << *expr << std::endl;
        visitor.init(inputs, *expr);

        auto f = visitor.apply(*expr);
        long double vars[3];

        // 相对于BCxx的偏移量
        // len_IC_start[idx]指示的是这个区域的所有位置的
        unsigned row_bais = len_IC_start[idx][regions[idx]->BCxx_loc[std::get<0>(loc)]];
        unsigned col_bais = len_IC_start[edge_idx][regions[edge_idx]->BCxx_loc[j]];

        for (int row = 0; row < row_num; row++)
        {
            for (int col = 0; col < col_num; col++)
            {
                long double h_val = static_cast<long double>(row + 1);
                long double n_val = static_cast<long double>(col + 1);

                BCxx.at(row + row_bais, col + col_bais) = integrate_dim0_param(f, vars, a, b, h_val, n_val);
                // BCxx[row + row_bais][col + col_bais] = integrate_dim0_param(f, vars, a, b, h_val, n_val);
            }
        }
    }
}

void AllRegions::splice_ES(Matrix2D &ES_xoz, Matrix2D &ES_yoz)
{
    unsigned row_s = 0;
    Matrix2D ES = Matrix2D(ES_xoz.rows() * ES_xoz.cols() + ES_yoz.rows() * ES_yoz.cols(), 1);

    for (int i = 0; i < all_region_num; i++)
    {
        const auto &idx_impl = regions[i]->impl;
        const auto &funcs = idx_impl->eq_ES;

        for (int j = 0; j < 6; j++)
        {
            // 如果当前这个取的边存在
            if (idx_impl->coeffs_exists[j])
            {
                RCP<const Basic> row_hn;
                RCP<const Basic> col_hn;
                unsigned row_num;
                unsigned col_num; // 行数和列数

                if (j <= 3)
                {
                    row_hn = idx_impl->h;
                    if (j == 0 || j == 2)
                        row_num = 1;
                    else
                        row_num = idx_impl->const_vals.H_max;
                }
                else
                {
                    row_hn = idx_impl->n;
                    row_num = idx_impl->const_vals.N_max;
                }

                // 判断是否是有源的
                if (funcs[j].size() != 0)
                {
                    Matrix2D ESx = Matrix2D(row_num, 1);
                    // 这个区域有源
                    for (auto &func : funcs[j])
                    {
                        Matrix2D ESxx = Matrix2D(row_num, 1);
                        auto expr = func.funcs;
                        // std::cout << *expr << std::endl;
                        double a = func.start_int;
                        double b = func.end_int;

                        // 转成std::function
                        SymEngine::vec_basic inputs;
                        inputs.push_back(func.x);
                        inputs.push_back(row_hn);

                        SymEngine::LambdaDoubleVisitor<long double> visitor;
                        // std::cout << *expr << std::endl;
                        visitor.init(inputs, *expr);

                        auto f = visitor.apply(*expr);
                        long double vars[3];

                        for (int row = 0; row < row_num; row++)
                        {
                            long double h_val = static_cast<long double>(row + 1);
                            long double n_val = static_cast<long double>(0);

                            ESxx.at(row, 0) = integrate_dim0_param(f, vars, a, b, h_val, n_val);
                        }
                        ESx += ESxx;
                    }

                    // 对应ESx找到对应的位置

                    ES.write_block(row_s, 0, ESx);
                    row_s += row_num;
                    // 处理row_sx
                }
                else
                {
                    // 区域无源，取0即可
                    row_s += row_num;
                }
            }
        }
    }
    ES.read_block(0, 0, ES_xoz.rows() * ES_xoz.cols(), 1, ES_xoz);
    ES.read_block(ES_xoz.rows() * ES_xoz.cols(), 0, ES.rows() * ES.cols(), 1, ES_yoz);
}

void AllRegions::split_IC(const std::vector<double> &IC_xoz, const std::vector<double> &IC_yoz)
{
}

void AllRegions::cal_IC_lens()
{
    this->all_len_IC.resize(this->region_num_xoz + this->region_num_yoz);
    this->len_IC_start.resize(this->region_num_xoz + this->region_num_yoz);

    for (int i = 0; i < this->region_num_xoz + this->region_num_yoz; i++)
    {
        const auto &idx_impl = this->regions[i]->impl;
        int len = 0;
        const auto &ce = idx_impl->coeffs_exists;
        int H = idx_impl->const_vals.H_max;
        int N = idx_impl->const_vals.N_max;
        if (ce[0])
        {
            len_IC_start[i].push_back(len);
            len += 1;
        }
        if (ce[1])
        {
            len_IC_start[i].push_back(len);
            len += H;
        }
        if (ce[2])
        {
            len_IC_start[i].push_back(len);
            len += 1;
        }
        if (ce[3])
        {
            len_IC_start[i].push_back(len);
            len += H;
        }
        if (ce[4])
        {
            len_IC_start[i].push_back(len);
            len += N;
        }
        if (ce[5])
        {
            len_IC_start[i].push_back(len);
            len += N;
        }

        this->all_len_IC[i] = len;
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

void AllRegions::input_calculate_area(double xl, double xr, double yb, double yt, double mu_r)
{
    this->divide_xoz->set_calculate_area(xl, xr, yb, yt, mu_r);
    this->divide_yoz->set_calculate_area(xl, xr, yb, yt, mu_r);
}

void AllRegions::pre_process()
{
    // 输入完所有的区域后，调用这个函数进行预处理
    this->divide_xoz->divide_regions();
    this->divide_xoz->findNeighbors();
    this->divide_xoz->cal_other_info();

    // 计算完成后，将数据输入到这个类中
    Plane_Info p1 = this->divide_xoz->rtn_Plane_Info(0);
    this->region_num_xoz = p1.regions_num;

    this->divide_yoz->divide_regions();
    this->divide_yoz->findNeighbors();
    this->divide_yoz->cal_other_info();

    Plane_Info p2 = this->divide_yoz->rtn_Plane_Info(this->region_num_xoz);
    this->region_num_yoz = p2.regions_num;

    // Merge xoz and yoz data
    this->all_region_num = this->region_num_xoz + this->region_num_yoz;

    auto append_regions = [&](const std::vector<Rect> *r1,
                              const std::vector<Rect> *r2)
    {
        regions_area.reserve(r1->size() + r2->size());

        regions_area.insert(regions_area.end(), r1->begin(), r1->end());
        regions_area.insert(regions_area.end(), r2->begin(), r2->end());
    };

    this->append_regions(this->regions_area, p1.divided_regions, p2.divided_regions);
    this->append_regions(this->current_regions, p1.current_regions, p2.current_regions);
    this->append_regions(this->all_casetypes, p1.all_casetype, p2.all_casetype);

    this->append_regions(this->all_H_max, p1.all_H_max, p2.all_H_max);
    this->append_regions(this->all_N_max, p1.all_N_max, p2.all_N_max);

    this->append_regions(this->all_mu_r, p1.all_mu_r, p2.all_mu_r);
    this->append_regions(this->all_J_r, p1.all_J_r, p2.all_J_r);

    this->append_regions(this->all_BC_types, p1.all_BC_types, p2.all_BC_types);

    this->append_regions(this->all_lefts, p1.all_lefts, p2.all_lefts);
    this->append_regions(this->all_rights, p1.all_rights, p2.all_rights);
    this->append_regions(this->all_tops, p1.all_tops, p2.all_tops);
    this->append_regions(this->all_bottoms, p1.all_bottoms, p2.all_bottoms);

    this->append_regions(this->current_regions_idx, p1.current_regions_idx, p2.current_regions_idx);
}

template <typename T>
void AllRegions::append_regions(std::vector<T> &out, const std::vector<T> *r1, const std::vector<T> *r2)
{
    out.reserve(r1->size() + r2->size());

    out.insert(out.end(), r1->begin(), r1->end());
    out.insert(out.end(), r2->begin(), r2->end());
}

template <typename Derived>
bool AllRegions::saveEigenToCSV(const std::string &filename,
                                const Eigen::MatrixBase<Derived> &mat,
                                char sep)
{
    std::ofstream file(filename);
    if (!file.is_open())
        return false;

    file << std::setprecision(17); // double / long double 安全精度

    for (int i = 0; i < mat.rows(); ++i)
    {
        for (int j = 0; j < mat.cols(); ++j)
        {
            file << mat(i, j);
            if (j + 1 < mat.cols())
                file << sep;
        }
        file << '\n';
    }

    return true;
}