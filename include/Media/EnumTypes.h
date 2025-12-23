#pragma once

#include <vector>
#include <array>
#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>

enum class CaseType
{
    AlleyAir,
    BTAir,
    NormalAir,
    FerriteCurrent,
    Aluminum
};

enum class BC_TYPE
{
    BBAA,
    AAAA,
    AABB
};

// 简单的矩形结构体
struct Rect
{
    float xl, xr, yb, yt;

    bool operator==(const Rect &other) const
    {
        return std::abs(xl - other.xl) < 1e-9 && std::abs(xr - other.xr) < 1e-9 &&
               std::abs(yb - other.yb) < 1e-9 && std::abs(yt - other.yt) < 1e-9;
    }
};

// 所有的这个plane的Region信息结构体
struct Plane_Info
{
    std::vector<Rect> *divided_regions;  // 划分好的区域
    int regions_num;                     // 区域数量
    std::vector<CaseType> *all_casetype; // 区域类型
    std::vector<Rect> *current_regions;  // 电流区域

    std::vector<int> *all_H_max;
    std::vector<int> *all_N_max;

    std::vector<float> *all_mu_r;
    std::vector<float> *all_J_r;

    // 边界邻接情况
    std::vector<std::vector<int>> *all_lefts;
    std::vector<std::vector<int>> *all_rights;
    std::vector<std::vector<int>> *all_tops;
    std::vector<std::vector<int>> *all_bottoms;

    std::vector<BC_TYPE> *all_BC_types;

    // 有源区域的下标
    std::vector<int> *current_regions_idx;
};

// 用于储存区域的常数变量
struct Region_Consts
{
    // 一共11个变量 2+5+4
    unsigned idx;
    unsigned all_regions_num; // 所有区域的数量

    Rect rect;
    unsigned H_max;
    unsigned N_max;
    float J_r;
    float mu_r;

    bool Ln;
    bool Rn;
    bool Tn;
    bool Bn;

    Region_Consts(const Region_Consts &c) : idx(c.idx), all_regions_num(c.all_regions_num), rect(c.rect),
                                            H_max(c.H_max), N_max(c.N_max),
                                            J_r(c.J_r), mu_r(c.mu_r),
                                            Ln(c.Ln), Rn(c.Rn), Tn(c.Tn), Bn(c.Bn)
    {
    }

    Region_Consts() = default;

    void set_Boundarys(bool Tn, bool Bn, bool Ln, bool Rn)
    {
        this->Tn = Tn;
        this->Bn = Bn;
        this->Ln = Ln;
        this->Rn = Rn;
    }
};

// 用于储存一个连接边界的所有方程
struct Boundary_Funcs
{
    std::array<SymEngine::RCP<const SymEngine::Basic>, 6> funcs; // 方程
    unsigned int idx;                                            // 对应邻接边的编号
    Rect rect;                                                   // 对应邻接区域
};

struct Integral_Func
{
    SymEngine::RCP<const SymEngine::Basic> funcs; // 积分函数
    SymEngine::RCP<const SymEngine::Basic> x;     // 积分变量

    float start_int; // 积分起始
    float end_int;   // 积分终止

    unsigned int edge_idx;

    Integral_Func(unsigned int edge_idx,
                  SymEngine::RCP<const SymEngine::Basic> &func,
                  SymEngine::RCP<const SymEngine::Basic> &x,
                  float start_int,
                  float end_int)
    {
        this->edge_idx = edge_idx;
        this->funcs = func;
        this->x = x;
        this->start_int = start_int;
        this->end_int = end_int;
    }

    Integral_Func() {}
};
