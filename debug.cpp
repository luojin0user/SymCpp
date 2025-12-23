#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/lambda_double.h>
#include <symengine/symbol.h>
#include <iostream>
#include <boost/math/quadrature/gauss.hpp>

using SymEngine::Expression;
using SymEngine::LambdaDoubleVisitor;
using SymEngine::vec_basic;

using boost::math::quadrature::gauss;

long double integrate_dim0(
    std::function<long double(const long double *)> &f,
    long double *x,
    long double a,
    long double b)
{
    boost::math::quadrature::gauss<long double, 15> quad;

    return quad.integrate(
        [&](long double t)
        {
            x[0] = t;
            return f(x);
        },
        a, b);
}

int main()
{
    // 1. 定义符号
    SymEngine::RCP<const SymEngine::Basic> x = SymEngine::symbol("x");
    SymEngine::RCP<const SymEngine::Basic> y = SymEngine::symbol("y");

    // 2. 构造表达式
    SymEngine::RCP<const SymEngine::Basic> expr = SymEngine::add(SymEngine::sinh(x), y);

    auto start_time = std::chrono::high_resolution_clock::now();

    // 3. 输入变量顺序（非常重要）
    vec_basic inputs;
    inputs.push_back(x);
    inputs.push_back(y);

    // 4. Visitor（long double）
    LambdaDoubleVisitor<long double> visitor;

    visitor.init(inputs, *expr);

    // 6. 直接生成 std::function
    auto f = visitor.apply(*expr);
    // f : std::function<long double(const long double*)>

    // 7. 输入值
    long double vals[] = {1000.0L};

    // 8. 执行
    long double result = f(vals);
    long double I = integrate_dim0(f, vals, 100, 100.1);

    // 6. 计时结束 (toc)
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Calculation finished in " << elapsed.count() << " seconds." << std::endl;

    std::cout << result << std::endl;
    std::cout << I;
}