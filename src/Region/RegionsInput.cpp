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

void RegionsInput::set_current_regions(float xl, float xr, float yb, float yt, float mu_r, float J_r)
{
    this->special_regions_area.push_back({xl, xr, yb, yt});
    this->special_regions_property.push_back({mu_r, J_r});
}

void RegionsInput::set_calculate_area(float xl, float xr, float yb, float yt, float mu_r)
{
    this->all_area = {xl, xr, yb, yt};
    this->air_mu_r = mu_r;
}

void RegionsInput::divide_regions()
{
    Rect S = this->all_area;
    this->special_regions_num = (int)this->special_regions_area.size();

    if (this->special_regions_area.empty())
    {
        this->divided_rects.clear();
        this->divided_rects.push_back(S);
        this->regions_num = 1;
        return;
    }

    // 1) 生成 y 边界（S 的上下 + 每个子区域的上下）
    std::vector<float> y_edges;
    y_edges.push_back(S.yb);
    y_edges.push_back(S.yt);
    for (const auto &sr : this->special_regions_area)
    {
        y_edges.push_back(sr.yb);
        y_edges.push_back(sr.yt);
    }

    // Unique and sort
    std::sort(y_edges.begin(), y_edges.end());
    auto last_y = std::unique(y_edges.begin(), y_edges.end(),
                              [](float a, float b)
                              { return std::abs(a - b) < 1e-10; });
    y_edges.erase(last_y, y_edges.end());

    // 2) 准备所有可能的 x 边界（用于中间需要分割的带）
    std::vector<float> all_x_edges;
    all_x_edges.push_back(S.xl);
    all_x_edges.push_back(S.xr);
    for (const auto &sr : this->special_regions_area)
    {
        all_x_edges.push_back(sr.xl);
        all_x_edges.push_back(sr.xr);
    }

    std::sort(all_x_edges.begin(), all_x_edges.end());
    auto last_x = std::unique(all_x_edges.begin(), all_x_edges.end(),
                              [](float a, float b)
                              { return std::abs(a - b) < 1e-10; });
    all_x_edges.erase(last_x, all_x_edges.end());

    std::vector<Rect> rects;

    // 3) 对每个 y 带处理
    for (size_t j = 0; j < y_edges.size() - 1; ++j)
    {
        float yB = y_edges[j];
        float yT = y_edges[j + 1];

        // 检查是否有子区域与该带交叠
        bool active = false;
        for (const auto &sr : this->special_regions_area)
        {
            // MATLAB: if (r(3) == yB) && (r(4) == yT)
            if (eq(sr.yb, yB) && eq(sr.yt, yT))
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
                float xL = all_x_edges[i];
                float xR = all_x_edges[i + 1];

                // 仅保留在 S 内的区段
                if (xL >= S.xl - 1e-10 && xR <= S.xr + 1e-10)
                {
                    // 当前区域如果是有源区域则需要特殊处理 (跳过，因为最后统一加)
                    bool is_special = false;
                    for (const auto &sr : this->special_regions_area)
                    {
                        if (eq(sr.xl, xL) && eq(sr.xr, xR) &&
                            eq(sr.yb, yB) && eq(sr.yt, yT))
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
    for (const auto &sr : this->special_regions_area)
    {
        rects.push_back(sr);
    }

    this->divided_rects = rects;
    this->regions_num = (int)rects.size();
}

bool RegionsInput::intervalsOverlap(const std::pair<float, float> &a, const std::pair<float, float> &b)
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
        float xiL = rects[i].xl;
        float xiR = rects[i].xr;
        float yiB = rects[i].yb;
        float yiT = rects[i].yt;

        for (int j = 0; j < N; ++j)
        {
            if (i == j)
                continue;

            float xjL = rects[j].xl;
            float xjR = rects[j].xr;
            float yjB = rects[j].yb;
            float yjT = rects[j].yt;

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
            this->all_BC_types[i] = BC_TYPE::BBAA;

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
            this->all_BC_types[i] = BC_TYPE::AABB;
            this->all_casetype[i] = CaseType::FerriteCurrent;
        }
        else
        {
            this->all_BC_types[i] = BC_TYPE::AAAA;
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
        this->all_mu_r_val[idx] = this->special_regions_property[i].first;
        this->all_J_r_val[idx] = this->special_regions_property[i].second;
    }
}

// --- Return functions ---

void RegionsInput::set_boundarys_idx(int start_idx)
{
    // 将src里面的所有数据都加上start_idx
    auto add_offset = [start_idx](std::vector<std::vector<int>> &src)
    {
        for (size_t i = 0; i < src.size(); i++)
        {
            for (int &val : src[i])
            {
                val += start_idx;
            }
        }
    };

    add_offset(this->all_lefts);
    add_offset(this->all_rights);
    add_offset(this->all_tops);
    add_offset(this->all_bottoms);
}

void RegionsInput::cal_current_idx(int start_idx)
{
    int i_start = this->regions_num - this->special_regions_num; // 0-based local index
    for (int i = 0; i < this->special_regions_num; ++i)
    {
        this->current_idx.push_back(i_start + i + start_idx);
    }
}

Plane_Info RegionsInput::rtn_Plane_Info(int start_idx)
{
    Plane_Info plain_info;

    plain_info.divided_regions = &this->divided_rects;
    plain_info.regions_num = this->regions_num;

    plain_info.current_regions = &this->special_regions_area;

    plain_info.all_H_max = &this->all_H_max;
    plain_info.all_N_max = &this->all_N_max;

    plain_info.all_mu_r = &this->all_mu_r_val;
    plain_info.all_J_r = &this->all_J_r_val;

    plain_info.all_BC_types = &this->all_BC_types;
    plain_info.all_casetype = &this->all_casetype;

    this->set_boundarys_idx(start_idx);
    plain_info.all_lefts = &this->all_lefts;
    plain_info.all_rights = &this->all_rights;
    plain_info.all_tops = &this->all_tops;
    plain_info.all_bottoms = &this->all_bottoms;

    this->cal_current_idx(start_idx);
    plain_info.current_regions_idx = &this->current_idx;

    return plain_info;
}