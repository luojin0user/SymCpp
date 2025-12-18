#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/lambda_double.h>
#include <symengine/symbol.h>
#include <iostream>

using SymEngine::Expression;
using SymEngine::LambdaDoubleVisitor;
using SymEngine::vec_basic;

int main()
{
    // 1. 定义符号
    Expression x("x");

    // 2. 构造表达式
    Expression expr = SymEngine::sinh(x) * x;

    // 3. 输入变量顺序（非常重要）
    vec_basic inputs;
    inputs.push_back(x.get_basic());

    // 4. Visitor（long double）
    LambdaDoubleVisitor<long double> visitor;

    visitor.init(inputs, *expr.get_basic());

    // 6. 直接生成 std::function
    auto f = visitor.apply(*expr.get_basic());
    // f : std::function<long double(const long double*)>

    // 7. 输入值
    long double vals[] = {1000.0L};

    // 8. 执行
    long double result = f(vals);

    std::cout << result << std::endl;
}
