#pragma once

// 引入必要的 SymEngine 头文件
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>
#include <vector>
#include <string>

// 为了不改变您的命名习惯，不使用额外的命名空间包裹类，直接定义
class BasicCase
{
public:
    // --- Properties (与 MATLAB 严格对应) ---
    SymEngine::RCP<const SymEngine::Symbol> x;
    SymEngine::RCP<const SymEngine::Symbol> y;

    int idx; // MATLAB: idx

    // --- 公共编号变量 ---
    SymEngine::RCP<const SymEngine::Symbol> xr;
    SymEngine::RCP<const SymEngine::Symbol> xl;
    SymEngine::RCP<const SymEngine::Symbol> yt;
    SymEngine::RCP<const SymEngine::Symbol> yl;

    // tau_x, tau_y 是计算结果，使用 Basic 类型
    SymEngine::RCP<const SymEngine::Basic> tau_x;
    SymEngine::RCP<const SymEngine::Basic> tau_y;

    // --- 公共函数 (作为成员变量存储的函数表达式) ---
    // 这里保留变量名，具体类型取决于它们在子类中是矩阵还是单个表达式
    // 暂时定义为 Basic，如果是矩阵后续可用 DenseMatrix
    // MATLAB: T_funcs, B_funcs...

    // --- 公共求和指标 ---
    SymEngine::RCP<const SymEngine::Symbol> h;
    SymEngine::RCP<const SymEngine::Symbol> n;

    // --- 公共参数 ---
    SymEngine::RCP<const SymEngine::Basic> beta_h;
    SymEngine::RCP<const SymEngine::Basic> lambda_n;
    SymEngine::RCP<const SymEngine::Basic> mu_0;
    SymEngine::RCP<const SymEngine::Symbol> mu_r; // 构造函数传入，通常视为符号

    int H_max;
    int N_max;

    // --- 边界积分项 ---
    SymEngine::RCP<const SymEngine::Symbol> c_0x;
    SymEngine::RCP<const SymEngine::Symbol> d_0x;
    SymEngine::RCP<const SymEngine::Basic> c_hx; // 可能会被置 0，所以用 Basic
    SymEngine::RCP<const SymEngine::Basic> d_hx;
    SymEngine::RCP<const SymEngine::Basic> e_ny;
    SymEngine::RCP<const SymEngine::Basic> f_ny;

    // 边界邻接区域逻辑变量
    bool Ln, Rn, Tn, Bn;

    int num_coeffs;
    // coeffs_exists 等数组逻辑建议在 C++ 中用 std::vector<int>

    // --- Methods ---
    BasicCase(int idx,
              const SymEngine::RCP<const SymEngine::Symbol> &xl,
              const SymEngine::RCP<const SymEngine::Symbol> &xr,
              const SymEngine::RCP<const SymEngine::Symbol> &yl,
              const SymEngine::RCP<const SymEngine::Symbol> &yt,
              int H_max, int N_max,
              const SymEngine::RCP<const SymEngine::Symbol> &mu_r);

    virtual ~BasicCase() = default;

    void apply_boundaries(bool Ln, bool Rn, bool Tn, bool Bn);

    // 虚函数，对应 MATLAB 中报错 "子类必须实现..."
    virtual void gen_solution_func(bool Ln, bool Rn, bool Tn, bool Bn);
    virtual void gen_coefficient_func();
};