#pragma once

#include <vector>
#include <string>
#include <map>

// SymEngine 基础头文件
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>

class RectangleCoil
{
public:
    // --- Properties (与 MATLAB 对应) ---

    // x和y方向的外长度
    double outer_len_x;
    double outer_len_y;

    // x和y方向的内长度
    double inner_len_x;
    double inner_len_y;

    // x,y,z中心位置，z是底面位置
    double x_loc;
    double y_loc;
    double z_loc;

    double mu0;
    int Nt;               // 总匝数
    int Nt_per_layer;     // 每一层的匝数
    int layer;            // 层数
    double wire_diameter; // 导线直径
    double thick;         // 厚度
    double Ir;            // 电流大小

    // 用于计算的区域
    double cal_xl;
    double cal_xr;
    double cal_yl;
    double cal_yr;
    double cal_zb;
    double cal_zt;

    // 修正因子 (SymEngine 表达式)
    SymEngine::RCP<const SymEngine::Basic> fx;
    SymEngine::RCP<const SymEngine::Basic> fy;

    // 分区属性
    double x_cl1, x_cl2, x_cr1, x_cr2;
    double m_xl, m_xr;
    double y_cl1, y_cl2, y_cr1, y_cr2;
    double m_yl, m_yr;
    double z_b, z_t;

    // 电流密度
    double Jr_xl, Jr_xr, Jr_yl, Jr_yr;

    // --- Methods ---

    // 构造函数
    RectangleCoil(double outer_len_x, double outer_len_y, int Nt_per_layer, int layer, double wire_diameter, double Ir);

    virtual ~RectangleCoil() = default;

    // 生成所有区域 (由于 AllRegions 类未定义，此处仅保留接口)
    void gen_all_regions();

    // 设置方形线圈的位置，x，y是中心位置，xoy坐标系下
    void set_Rcoil_loc(double x, double y, double z);

    // 设置计算区域并触发因子计算
    void set_calculate_area(double xl, double xr, double yl, double yr, double zb, double zt, double mu0);

    // 用于2D到3D的修正因子计算 (x方向)
    void factor_x();

    // 辅助函数：计算 factor_x 的具体数值
    double factor_x_real(double x_val, double y_val, double z_val);

    // 用于2D到3D的修正因子计算 (y方向)
    void factor_y();

    // 绘图函数 plot3D 和 drawBox 已根据指示移除
};