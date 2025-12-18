#pragma once

#include "BasicCase.hpp"
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>
#include <vector>
#include <memory>
#include <map>

class Case1 : public BasicCase
{
public:
    // --- 补充声明 MATLAB 代码中使用的成员变量 ---
    // 这些变量在提供的 BasicCase.m 中未显式定义，但在 Case1 中被大量使用
    // 因此在这里定义以保证逻辑完整性。

    // 区域连接关系
    std::vector<int> tops;
    std::vector<int> bottoms;
    std::vector<int> lefts;
    std::vector<int> rights;

    // 所有区域的指针 (模拟 obj.all_regions)
    std::vector<std::shared_ptr<BasicCase>> all_regions;

    // 边界函数输入 (模拟 B_funcs, T_funcs 等 cell 数组)
    std::vector<SymEngine::RCP<const SymEngine::Basic>> B_funcs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> T_funcs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> L_funcs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> R_funcs;

    // 有源项函数 (ES_funcs)
    std::vector<SymEngine::RCP<const SymEngine::Basic>> T_ESfuncs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> B_ESfuncs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> L_ESfuncs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> R_ESfuncs;

    // 有源区域标记
    std::vector<bool> ES_regions;

    // 边界方程存储 (模拟 cell 数组，使用嵌套 vector 存储方程)
    // 结构: [row][col] -> Equation
    std::vector<std::vector<SymEngine::RCP<const SymEngine::Basic>>> eq_c_hx;
    std::vector<std::vector<SymEngine::RCP<const SymEngine::Basic>>> eq_d_hx;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_e_ny; // 1D cell
    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_f_ny; // 1D cell
    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_ES;

    // 系数存在性 (coeffs_exists)
    std::vector<int> coeffs_exists;

    // 位置映射 (BCfuncs_loc_map) - 假设大小固定或足够大
    // Key: column index (region idx), Value: pair(type, func_num)
    std::map<int, std::pair<int, int>> BCfuncs_loc_map;

    // 磁势和磁场表达式
    SymEngine::RCP<const SymEngine::Basic> A_zx_expr;
    SymEngine::RCP<const SymEngine::Basic> A_zy_expr;

    SymEngine::RCP<const SymEngine::Basic> B_x_x;
    SymEngine::RCP<const SymEngine::Basic> B_x_y;
    SymEngine::RCP<const SymEngine::Basic> B_y_x;
    SymEngine::RCP<const SymEngine::Basic> B_y_y;

    // 最终方程
    SymEngine::RCP<const SymEngine::Basic> eq_A_z;
    SymEngine::RCP<const SymEngine::Basic> eq_B_x;
    // -------------------------------------------------------

    // 构造函数
    Case1(int idx,
          const SymEngine::RCP<const SymEngine::Symbol> &xl,
          const SymEngine::RCP<const SymEngine::Symbol> &xr,
          const SymEngine::RCP<const SymEngine::Symbol> &yl,
          const SymEngine::RCP<const SymEngine::Symbol> &yt,
          bool Ln, bool Rn, bool Tn, bool Bn,
          int H_max, int N_max,
          const SymEngine::RCP<const SymEngine::Symbol> &mu_r);

    virtual ~Case1() = default;

    // 重写父类虚函数
    void gen_solution_func(bool Ln, bool Rn, bool Tn, bool Bn) override; // 注意：父类定义中包含参数，虽此处MATLAB代码未显式用
    void gen_solution_func();                                            // MATLAB 代码中的无参版本

    void gen_coefficient_func() override;
};