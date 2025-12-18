#include <iostream>
#include <chrono>
#include <cmath>

// 引入之前生成的头文件
#include "RectangleCoil.hpp"
#include "AllRegions.hpp"

// 定义数学常数 (若未定义)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main()
{
    // 1. 计时开始 (tic)
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "Starting simulation..." << std::endl;

    // 2. 初始化 RectangleCoil
    // MATLAB: coil = RectangleCoil(0.08, 0.18, 28, 56, 0.707e-3, 5);
    double outer_len_x = 0.08;
    double outer_len_y = 0.18;
    int Nt_per_layer = 28;
    int layer = 56;
    double wire_diameter = 0.707e-3;
    double Ir = 5.0;

    RectangleCoil coil(outer_len_x, outer_len_y, Nt_per_layer, layer, wire_diameter, Ir);

    // 3. 设置线圈位置
    // MATLAB: coil.set_Rcoil_loc(0.14, 0.14, 0.03);
    coil.set_Rcoil_loc(0.14, 0.14, 0.03);

    // 4. 设置计算区域
    // MATLAB: coil.set_calculate_area(0, 0.28, 0, 0.28, 0, 0.24, 4*pi*1e-7);
    double mu0 = 4 * M_PI * 1e-7;
    coil.set_calculate_area(0.0, 0.28, 0.0, 0.28, 0.0, 0.24, mu0);

    coil.gen_all_regions();

    // 6. 计时结束 (toc)
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Calculation finished in " << elapsed.count() << " seconds." << std::endl;

    return 0;
}