#pragma once

#include "BasicCase.hpp"
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>
#include <vector>
#include <memory>
#include <map>

class Case2 : public BasicCase
{
public:
    // --- Properties 特有属性 ---
    SymEngine::RCP<const SymEngine::Basic> J_z;
    SymEngine::RCP<const SymEngine::Basic> A_y_P;
    SymEngine::RCP<const SymEngine::Basic> B_x_P;
    SymEngine::RCP<const SymEngine::Basic> B_y_P;

    // 方程对象
    SymEngine::RCP<const SymEngine::Basic> eq_A_y_P;
    SymEngine::RCP<const SymEngine::Basic> eq_B_x_P;
    SymEngine::RCP<const SymEngine::Basic> eq_B_y_P;

    // --- 补充声明 MATLAB 代码中使用的成员变量 (模拟环境) ---
    // 假设这些在运行时由外部填充
    std::vector<int> tops;
    std::vector<int> bottoms;
    std::vector<int> lefts;
    std::vector<int> rights;

    // 所有区域 (模拟 obj.all_regions)
    std::vector<std::shared_ptr<BasicCase>> all_regions;

    // 边界函数输入
    std::vector<SymEngine::RCP<const SymEngine::Basic>> B_funcs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> T_funcs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> L_funcs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> R_funcs;

    // 有源项函数
    std::vector<SymEngine::RCP<const SymEngine::Basic>> T_ESfuncs;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> B_ESfuncs;

    // 有源区域标记
    std::vector<bool> ES_regions;

    // 系数存在性
    std::vector<int> coeffs_exists;

    // 方程存储容器
    // 注意：Case2 中 eq_c_hx 在 MATLAB 代码中是一维索引 eq_c_hx{i}，与 Case1 不同
    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_c_hx;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_d_hx;

    // 直流分量/常数项
    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_c0x;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_d0x;

    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_e_ny;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_f_ny;
    std::vector<SymEngine::RCP<const SymEngine::Basic>> eq_ES;

    // 位置映射
    std::map<int, std::pair<int, int>> BCfuncs_loc_map;

    // 求解结果暂存
    SymEngine::RCP<const SymEngine::Basic> A_zx_expr;
    SymEngine::RCP<const SymEngine::Basic> A_zy_expr;
    SymEngine::RCP<const SymEngine::Basic> B_x_x;
    SymEngine::RCP<const SymEngine::Basic> B_x_y;
    SymEngine::RCP<const SymEngine::Basic> B_y_x;
    SymEngine::RCP<const SymEngine::Basic> B_y_y;

    SymEngine::RCP<const SymEngine::Basic> eq_A_z;
    SymEngine::RCP<const SymEngine::Basic> eq_B_x;

    // --- Methods ---
    Case2(int idx,
          const SymEngine::RCP<const SymEngine::Symbol> &xl,
          const SymEngine::RCP<const SymEngine::Symbol> &xr,
          const SymEngine::RCP<const SymEngine::Symbol> &yl,
          const SymEngine::RCP<const SymEngine::Symbol> &yt,
          bool Ln, bool Rn, bool Tn, bool Bn,
          int H_max, int N_max,
          const SymEngine::RCP<const SymEngine::Symbol> &mu_r,
          const SymEngine::RCP<const SymEngine::Basic> &J_r); // J_r 可以是数值或符号

    virtual ~Case2() = default;

    void gen_solution_func(bool Ln, bool Rn, bool Tn, bool Bn) override; // 适配父类虚函数

    void gen_coefficient_func() override;
};