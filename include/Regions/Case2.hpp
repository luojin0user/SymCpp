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

class Case2 : public BasicCase
{
public:
    // --- Properties 特有属性 ---
    SymEngine::RCP<const SymEngine::Basic> J_z;
    SymEngine::RCP<const SymEngine::Basic> A_y_P;
    SymEngine::RCP<const SymEngine::Basic> B_x_P;
    SymEngine::RCP<const SymEngine::Basic> B_y_P;

    // 构造函数
    Case2(Region_Consts &c);

    virtual ~Case2() = default;

    // 重写父类虚函数
    void gen_solution_func();

    void gen_coefficient_func();
};