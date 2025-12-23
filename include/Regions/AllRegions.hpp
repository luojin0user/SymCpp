#pragma once

#include "EnumTypes.h"
#include "RegionsInput.hpp"
#include "Matrix2D.hpp"

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <algorithm>

// SymEngine
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

#include <boost/math/quadrature/gauss.hpp>

#include <Eigen/Dense>

// Forward Declarations
class Region;
class RegionsInput; // Assuming this class exists based on MATLAB code
enum class CaseType;
enum class BC_TYPE;
struct Rect;

class AllRegions
{
public:
    // --- Properties ---

    // 存储所有区域 (Use shared_ptr for automatic memory management)
    // std::vector<std::shared_ptr<Region>> regions;
    std::vector<std::unique_ptr<Region>> regions;

    // 用于处理输入区域等 (Assuming RegionsInput is a defined class)
    std::unique_ptr<RegionsInput> divide_xoz;
    std::unique_ptr<RegionsInput> divide_yoz;

    int region_num_xoz; // 区域数量
    int region_num_yoz; // 区域数量
    int all_region_num;

    std::vector<int> all_H_max;
    std::vector<int> all_N_max;

    std::vector<float> all_mu_r;
    std::vector<float> all_J_r;

    std::vector<int> all_len_IC;
    std::vector<std::vector<unsigned>> len_IC_start;

    // 区域，每一行是一个矩形区域 [xl, xr, yl, yt]
    std::vector<Rect> regions_area;

    std::vector<BC_TYPE> all_BC_types;
    std::vector<CaseType> all_casetypes;

    // 有源区域
    std::vector<Rect> current_regions;
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

    void cal_all_BCs();

    void splice_BC(Matrix2D &BC_xoz, Matrix2D &BC_yoz);

    void splice_BCxx(int idx, int edge_idx, Matrix2D &BCxx);

    // 拼接ES矩阵
    void splice_ES(Matrix2D &ES_xoz, Matrix2D &ES_yoz);

    void split_IC(const std::vector<double> &IC_xoz, const std::vector<double> &IC_yoz);

    void cal_IC_lens();

    // Calculation Helpers
    // void cal_Bx_By_region(const std::string &plane, const SplitICResult &ICs, int region_num,
    //                       const std::vector<double> &x0, const std::vector<double> &y0,
    //                       std::vector<double> &Bx, std::vector<double> &By);

    // void cal_Bx_By(const std::string &plane, const SplitICResult &ICs,
    //                const std::vector<double> &x0, const std::vector<double> &y0,
    //                std::vector<double> &Bx_out, std::vector<double> &By_out);

    // Input proxies
    void input_current_region(const std::string &plain, double xl, double xr, double yb, double yt, double mu_r, double J_r);
    void input_calculate_area(double xl, double xr, double yb, double yt, double mu_r);

    void pre_process();

    template <typename T>
    void append_regions(std::vector<T> &out, const std::vector<T> *r1, const std::vector<T> *r2);

    template <typename Derived>
    bool saveEigenToCSV(const std::string &filename,
                        const Eigen::MatrixBase<Derived> &mat,
                        char sep = ',');

private:
    inline long double integrate_dim0_param(
        const std::function<long double(const long double *)> &f,
        long double *x,
        long double a,
        long double b,
        long double h,
        long double n)
    {
        boost::math::quadrature::gauss<long double, 15> quad;

        return quad.integrate(
            [&](long double t)
            {
                x[0] = t; // x
                x[1] = h; // h
                x[2] = n; // n
                return f(x);
            },
            a, b);
    }
};