#include <iostream>
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/printers.h>
#include <symengine/basic.h>

// 使用 SymEngine 的命名空间
using SymEngine::add;
using SymEngine::Basic;
using SymEngine::expand;
using SymEngine::integer;
using SymEngine::mul;
using SymEngine::pow;
using SymEngine::RCP;
using SymEngine::symbol;

#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <iostream>

using SymEngine::Expression;
using SymEngine::symbol;

int main()
{
    // 定义符号
    Expression x("x");
    Expression y("y");

    // 构造表达式
    Expression expr = x * x + 2 * x * y + y * y;

    // 输出表达式
    std::cout << "expr = " << expr << std::endl;

    // 求导
    Expression dex_dx = expr.diff(x);
    Expression dex_dy = expr.diff(y);

    std::cout << "d(expr)/dx = " << dex_dx << std::endl;
    std::cout << "d(expr)/dy = " << dex_dy << std::endl;

    // 代值计算
    SymEngine::map_basic_basic subs_map;
    subs_map[x] = integer(1);
    subs_map[y] = integer(2);

    Expression value = expr.subs(subs_map);
    std::cout << "expr(x=1,y=2) = " << value << std::endl;

    return 0;
}
