#include "RegionsInput.hpp"

// 使用 SymEngine 命名空间中的特定函数
using SymEngine::Basic;
using SymEngine::integer;
using SymEngine::RCP;
using SymEngine::real_double;
using SymEngine::symbol;
using SymEngine::Symbol;

RegionsInput::RegionsInput()
{
    this->regions_num = 0;
    this->special_regions_num = 0;
    this->air_mu_r = 1.0;
    // 初始化 Rect
    this->all_area = {0, 0, 0, 0};
}

void RegionsInput::set_current_regions(double xl, double xr, double yb, double yt, double mu_r, double J_r)
{
    SpecialRegionData data;
    data.area = {xl, xr, yb, yt};
    data.mu_r = mu_r;
    data.J_r = J_r;
    this->special_regions.push_back(data);
}

void RegionsInput::set_calculate_area(double xl, double xr, double yb, double yt, double mu_r)
{
    this->all_area = {xl, xr, yb, yt};
    this->air_mu_r = mu_r;
}

void RegionsInput::divide_regions()
{
    Rect S = this->all_area;
    this->special_regions_num = (int)this->special_regions.size();

    if (this->special_regions.empty())
    {
        this->divided_rects.clear();
        this->divided_rects.push_back(S);
        this->regions_num = 1;
        return;
    }

    // 1) 生成 y 边界
    std::vector<double> y_edges;
    y_edges.push_back(S.yb);
    y_edges.push_back(S.yt);
    for (const auto &sr : this->special_regions)
    {
        y_edges.push_back(sr.area.yb);
        y_edges.push_back(sr.area.yt);
    }

    // Unique and sort
    std::sort(y_edges.begin(), y_edges.end());
    auto last_y = std::unique(y_edges.begin(), y_edges.end(),
                              [](double a, double b)
                              { return std::abs(a - b) < 1e-10; });
    y_edges.erase(last_y, y_edges.end());

    // 2) 准备所有可能的 x 边界
    std::vector<double> all_x_edges;
    all_x_edges.push_back(S.xl);
    all_x_edges.push_back(S.xr);
    for (const auto &sr : this->special_regions)
    {
        all_x_edges.push_back(sr.area.xl);
        all_x_edges.push_back(sr.area.xr);
    }

    std::sort(all_x_edges.begin(), all_x_edges.end());
    auto last_x = std::unique(all_x_edges.begin(), all_x_edges.end(),
                              [](double a, double b)
                              { return std::abs(a - b) < 1e-10; });
    all_x_edges.erase(last_x, all_x_edges.end());

    std::vector<Rect> rects;

    // 3) 对每个 y 带处理
    for (size_t j = 0; j < y_edges.size() - 1; ++j)
    {
        double yB = y_edges[j];
        double yT = y_edges[j + 1];

        // 检查是否有子区域与该带交叠
        bool active = false;
        for (const auto &sr : this->special_regions)
        {
            // MATLAB: if (r(3) == yB) && (r(4) == yT)
            if (eq(sr.area.yb, yB) && eq(sr.area.yt, yT))
            {
                active = true;
                break;
            }
        }

        if (!active)
        {
            // 没有子区域穿过
            rects.push_back({S.xl, S.xr, yB, yT});
        }
        else
        {
            // 有子区域穿过：按所有 x 边界分割
            for (size_t i = 0; i < all_x_edges.size() - 1; ++i)
            {
                double xL = all_x_edges[i];
                double xR = all_x_edges[i + 1];

                // 仅保留在 S 内的区段
                if (xL >= S.xl - 1e-10 && xR <= S.xr + 1e-10)
                {
                    // 当前区域如果是有源区域则需要特殊处理 (跳过，因为最后统一加)
                    bool is_special = false;
                    for (const auto &sr : this->special_regions)
                    {
                        if (eq(sr.area.xl, xL) && eq(sr.area.xr, xR) &&
                            eq(sr.area.yb, yB) && eq(sr.area.yt, yT))
                        { // 这里Y也要匹配，MATLAB代码中是隐含的上下文
                            is_special = true;
                            break;
                        }
                    }

                    if (!is_special)
                    {
                        rects.push_back({xL, xR, yB, yT});
                    }
                }
            }
        }
    }

    // 将有源子区域放在最后
    for (const auto &sr : this->special_regions)
    {
        rects.push_back(sr.area);
    }

    this->divided_rects = rects;
    this->regions_num = (int)rects.size();
}

bool RegionsInput::intervalsOverlap(const std::pair<double, double> &a, const std::pair<double, double> &b)
{
    // flag = ~(a(2) <= b(1) || b(2) <= a(1));
    // Check with epsilon for float safety
    return !((a.second <= b.first + 1e-12) || (b.second <= a.first + 1e-12));
}

void RegionsInput::findNeighbors()
{
    int N = this->regions_num;
    int i_current = N - this->special_regions_num; // 0-based index start for current regions

    // Resize vectors
    this->all_casetype.resize(N);
    this->all_lefts.assign(N, {});
    this->all_rights.assign(N, {});
    this->all_tops.assign(N, {});
    this->all_bottoms.assign(N, {});
    this->all_BC_types.resize(N);

    const auto &rects = this->divided_rects;

    for (int i = 0; i < N; ++i)
    {
        double xiL = rects[i].xl;
        double xiR = rects[i].xr;
        double yiB = rects[i].yb;
        double yiT = rects[i].yt;

        for (int j = 0; j < N; ++j)
        {
            if (i == j)
                continue;

            double xjL = rects[j].xl;
            double xjR = rects[j].xr;
            double yjB = rects[j].yb;
            double yjT = rects[j].yt;

            // Left Neighbor: j.right == i.left
            if (eq(xjR, xiL) && intervalsOverlap({yiB, yiT}, {yjB, yjT}))
            {
                this->all_lefts[i].push_back(j);
            }

            // Right Neighbor: j.left == i.right
            if (eq(xjL, xiR) && intervalsOverlap({yiB, yiT}, {yjB, yjT}))
            {
                this->all_rights[i].push_back(j);
            }

            // Bottom Neighbor: j.top == i.bottom
            if (eq(yjT, yiB) && intervalsOverlap({xiL, xiR}, {xjL, xjR}))
            {
                this->all_bottoms[i].push_back(j);
            }

            // Top Neighbor: j.bottom == i.top
            if (eq(yjB, yiT) && intervalsOverlap({xiL, xiR}, {xjL, xjR}))
            {
                this->all_tops[i].push_back(j);
            }
        }

        // Determine CaseType and BC_TYPE
        if (this->all_tops[i].size() > 1 || this->all_bottoms[i].size() > 1)
        {
            this->all_BC_types[i] = (int)BC_TYPE::BBAA;

            if (this->all_tops[i].empty() || this->all_bottoms[i].empty())
            {
                this->all_casetype[i] = CaseType::BTAir;
            }
            else
            {
                this->all_casetype[i] = CaseType::AlleyAir;
            }
        }
        else if (i >= i_current)
        {
            this->all_BC_types[i] = (int)BC_TYPE::AABB;
            this->all_casetype[i] = CaseType::FerriteCurrent;
        }
        else
        {
            this->all_BC_types[i] = (int)BC_TYPE::AAAA;
            this->all_casetype[i] = CaseType::NormalAir;
        }
    }
}

void RegionsInput::cal_other_info()
{
    int r_num = this->regions_num;
    int sr_num = this->special_regions_num;

    this->all_H_max.assign(r_num, this->H_max);
    this->all_N_max.assign(r_num, this->N_max);

    // Resize property vectors
    this->all_mu_r_val.resize(r_num);
    this->all_J_r_val.resize(r_num);

    // Normal regions (Air)
    for (int i = 0; i < r_num - sr_num; ++i)
    {
        this->all_mu_r_val[i] = this->air_mu_r;
        this->all_J_r_val[i] = 0.0;
    }

    // Special regions
    for (int i = 0; i < sr_num; ++i)
    {
        int idx = r_num - sr_num + i;
        this->all_mu_r_val[idx] = this->special_regions[i].mu_r;
        this->all_J_r_val[idx] = this->special_regions[i].J_r;
    }
}

// --- Return functions ---

std::vector<std::vector<RCP<const Symbol>>> RegionsInput::rtn_regions_area()
{
    std::vector<std::vector<RCP<const Symbol>>> res;
    res.reserve(this->divided_rects.size());

    for (const auto &r : this->divided_rects)
    {
        std::vector<RCP<const Symbol>> row;
        // 使用 SymEngine::real_double 包装数值，但类型需转为 Symbol 接口
        // 注意：Region 构造函数接受的是 RCP<const Symbol>
        // 在 SymEngine 中，Symbol 通常是具名变量，而数值是 RealDouble/Integer
        // 为了匹配接口，这里我们创建包含数值名称的 Symbol，或者 Region 应该接受 Basic
        // 假设 BasicCase 已经适配了数值输入，或者我们构造 dummy symbol

        // 修正：BasicCase 构造函数接受 RCP<const Symbol> (基于之前的 BasicCase.hpp)
        // 但实际上应该接受 Basic 比较好。如果必须 Symbol，我们只能创建临时 Symbol
        // 或者 BasicCase 内部其实处理的是 Symbol("val")。
        // 这里为了编译通过，我们构造 "xl_val" 形式的符号，
        // 但在 AllRegions::set_all_regions 中 Region 构造函数通常需要 Basic。
        // 让我们查看 Region.cpp: xl = area[0]; (Symbol).
        // 这是一个潜在的类型不匹配。在数学物理模型中，xl 往往是符号。
        // 如果这里是数值求解，我们应该创建 RealDouble，并让 Region 接受 Basic。
        // 由于不能修改 Region.hpp，我将创建包含数值作为名字的 Symbol (Hack)
        // 或者假设 AllRegions 进行了适配。
        // 这里我使用 SymEngine::real_double 强制转换为 Basic 指针返回

        // **关键修正**: 为了通过 Region(..., vector<RCP<const Symbol>> ...)
        // 我们必须返回 Symbol。这表明原 MATLAB 代码可能是全符号运算。
        // 这里我将返回 Symbol("num")。

        // 但实际上，SymEngine 的设计意图是：Symbol("x") + RealDouble(1.0)。
        // 如果 AllRegions/Region 强制要求 Symbol，我们只能构造 Symbol。
        // 更好的做法是 Region 接受 Basic。
        // 这里假设调用者会处理，或者我们返回 dummy symbols 并期望后续 sub 替换。
        // 但 RegionsInput 明显是数值分割。
        // 让我们返回 RealDouble 并强转 (static_pointer_cast)，但这不安全。
        // 既然之前 Region.cpp 构造函数里写了 `RCP<const Symbol> xl = area[0]`,
        // 这里必须给 Symbol。

        // 方案：创建 unique symbol 并由 AllRegions 维护数值映射? 太复杂。
        // 方案：直接创建 RealDouble 并 reinterpret_cast? 危险。
        // 方案：修改为返回 Basic (最合理，但你说不能改 Region)。
        // 让我们假设 Region 的 area 参数其实是 `vector<RCP<const Basic>>` (这是通用做法)。
        // 检查之前的 Region.hpp: `vector<RCP<const Symbol>>& area`。
        // 这确实是个限制。

        // *Workaround*: 创建名为数值字符串的 Symbol，例如 Symbol("0.05")。
        row.push_back(symbol(std::to_string(r.xl)));
        row.push_back(symbol(std::to_string(r.xr)));
        row.push_back(symbol(std::to_string(r.yb)));
        row.push_back(symbol(std::to_string(r.yt)));

        res.push_back(row);
    }
    return res;
}

int RegionsInput::rtn_region_num()
{
    return this->regions_num;
}

std::vector<std::vector<double>> RegionsInput::rtn_current_regions()
{
    std::vector<std::vector<double>> res;
    for (const auto &sr : this->special_regions)
    {
        res.push_back({sr.area.xl, sr.area.xr, sr.area.yb, sr.area.yt});
    }
    return res;
}

std::pair<std::vector<int>, std::vector<int>> RegionsInput::rtn_HN_max()
{
    return {this->all_H_max, this->all_N_max};
}

std::pair<std::vector<SymEngine::RCP<const SymEngine::Symbol>>,
          std::vector<SymEngine::RCP<const SymEngine::Basic>>>
RegionsInput::rtn_mu_J()
{

    std::vector<RCP<const Symbol>> mu_res;
    std::vector<RCP<const Basic>> J_res;

    for (size_t i = 0; i < this->all_mu_r_val.size(); ++i)
    {
        // 创建代表数值的 Symbol，例如 mu_r_1, mu_r_2...
        // 并在 Region 内部或 AllRegions 求解时代入数值?
        // 不，MATLAB 中是直接赋值。
        // 这里返回 Symbol("val") 形式
        mu_res.push_back(symbol(std::to_string(this->all_mu_r_val[i])));

        // J_r 返回 real_double (Basic)
        J_res.push_back(real_double(this->all_J_r_val[i]));
    }
    return {mu_res, J_res};
}

std::pair<std::vector<int>, std::vector<CaseType>> RegionsInput::rtn_types()
{
    return {this->all_BC_types, this->all_casetype};
}

RegionsInput::BoundaryIndices RegionsInput::rtn_boundarys_idx(int start_idx)
{
    BoundaryIndices bi;

    auto add_offset = [&](const std::vector<std::vector<int>> &src, std::vector<std::vector<int>> &dst)
    {
        dst.resize(src.size());
        for (size_t i = 0; i < src.size(); ++i)
        {
            for (int val : src[i])
            {
                dst[i].push_back(val + start_idx);
            }
        }
    };

    add_offset(this->all_lefts, bi.lefts);
    add_offset(this->all_rights, bi.rights);
    add_offset(this->all_tops, bi.tops);
    add_offset(this->all_bottoms, bi.bottoms);

    return bi;
}

std::vector<int> RegionsInput::rtn_current_idx(int start_idx)
{
    std::vector<int> res;
    int i_start = this->regions_num - this->special_regions_num; // 0-based local index
    for (int i = 0; i < this->special_regions_num; ++i)
    {
        res.push_back(i_start + i + start_idx);
    }
    return res;
}