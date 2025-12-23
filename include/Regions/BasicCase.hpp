#pragma once

// 引入必要的 SymEngine 头文件
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>
#include <vector>
#include <string>
#include <array>

#include "EnumTypes.h"

class BasicCase
{
public:
    const Region_Consts const_vals;

    // --- Properties (与 MATLAB 严格对应) ---
    SymEngine::RCP<const SymEngine::Symbol> x;
    SymEngine::RCP<const SymEngine::Symbol> y;

    int idx; // MATLAB: idx

    // --- 公共编号变量 ---
    SymEngine::RCP<const SymEngine::Basic> xr;
    SymEngine::RCP<const SymEngine::Basic> xl;
    SymEngine::RCP<const SymEngine::Basic> yt;
    SymEngine::RCP<const SymEngine::Basic> yl;

    // tau_x, tau_y 是计算结果，使用 Basic 类型
    SymEngine::RCP<const SymEngine::Basic> tau_x;
    SymEngine::RCP<const SymEngine::Basic> tau_y;

    // --- 公共求和指标 ---
    SymEngine::RCP<const SymEngine::Basic> h;
    SymEngine::RCP<const SymEngine::Basic> n;

    // --- 公共参数 ---
    SymEngine::RCP<const SymEngine::Basic> beta_h;
    SymEngine::RCP<const SymEngine::Basic> lambda_n;
    SymEngine::RCP<const SymEngine::Basic> mu_0;
    SymEngine::RCP<const SymEngine::Basic> mu_r; // 构造函数传入，通常视为符号

    // --- 边界积分项 ---
    SymEngine::RCP<const SymEngine::Basic> c_0x;
    SymEngine::RCP<const SymEngine::Basic> d_0x;
    SymEngine::RCP<const SymEngine::Basic> c_hx;
    SymEngine::RCP<const SymEngine::Basic> d_hx;
    SymEngine::RCP<const SymEngine::Basic> e_ny;
    SymEngine::RCP<const SymEngine::Basic> f_ny;

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

    // 有源区域标记，有源为true
    std::vector<bool> ES_regions;

    int num_coeffs;
    // 系数存在性 (coeffs_exists)
    std::vector<int> coeffs_exists;

    // 磁势和磁场表达式
    SymEngine::RCP<const SymEngine::Basic> A_zx_expr;
    SymEngine::RCP<const SymEngine::Basic> A_zy_expr;

    SymEngine::RCP<const SymEngine::Basic> B_x_x;
    SymEngine::RCP<const SymEngine::Basic> B_x_y;
    SymEngine::RCP<const SymEngine::Basic> B_y_x;
    SymEngine::RCP<const SymEngine::Basic> B_y_y;

    /***********下面这两个是计算的结果，会被外部调用*********** */

    // 边界方程存储 (使用嵌套 vector 存储方程)
    // 结构: [row][col] -> Equation
    std::array<std::vector<std::array<Integral_Func, 6>>, 6> eq_BC;

    // 用于储存邻接区域edge_idx对应上面方程的坐标位置，eq_BC_loc下标即为edge_idx，值为坐标，bool代表是否有c0/d0
    std::vector<std::tuple<int, int, bool>> eq_BC_loc;

    // 所有的ES，可能有6个边界（当前区域的c0 c d0 d e f），内部可能有多个，需要后面进行相加，外面则是根据边界分开的
    std::array<std::vector<Integral_Func>, 6> eq_ES;

    // --- Methods ---
    BasicCase(Region_Consts &c);

    virtual ~BasicCase() = default;

    // 虚函数
    virtual void gen_solution_func() = 0;
    virtual void gen_coefficient_func() = 0;

    void gen_integral(const std::vector<Boundary_Funcs> &BC_func,
                      const std::vector<SymEngine::RCP<const SymEngine::Basic>> &BC_ESfunc,
                      int BTLR, // 指示上下左右，1开始
                      std::vector<std::array<Integral_Func, 6UL>> &eq,
                      std::vector<Integral_Func> &eq_ES,
                      std::vector<std::array<Integral_Func, 6UL>> *eq_c0d0 = nullptr);

private:
    void apply_boundaries();
};