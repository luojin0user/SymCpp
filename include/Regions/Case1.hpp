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
};