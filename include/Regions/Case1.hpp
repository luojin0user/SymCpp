#pragma once

#include "BasicCase.hpp"
#include "Region.hpp"
#include "EnumTypes.h"

#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>
#include <vector>
#include <memory>
#include <map>
#include <array>

class Case1 : public BasicCase
{
public:
    // 最终方程
    SymEngine::RCP<const SymEngine::Basic> eq_A_z;
    SymEngine::RCP<const SymEngine::Basic> eq_B_x;
    // -------------------------------------------------------

    // 构造函数
    Case1(Region_Consts &c);

    virtual ~Case1() = default;

    // 重写父类虚函数
    void gen_solution_func();

    void gen_coefficient_func();

    void gen_integral(const std::vector<Boundary_Funcs> &BC_func,
                      const std::vector<SymEngine::RCP<const SymEngine::Basic>> &BC_ESfunc,
                      int BTLR, // 指示上下左右，1开始
                      std::vector<std::array<Integral_Func, 6UL>> &eq,
                      std::vector<Integral_Func> &eq_ES,
                      std::vector<std::array<Integral_Func, 6UL>> *eq_c0d0 = nullptr);
};