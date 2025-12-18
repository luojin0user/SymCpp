#pragma once

#include <vector>
#include <cmath>
#include <iostream>

// SymEngine
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/real_double.h>
#include <symengine/eval_double.h>

namespace SplitCurveUtil
{

    using SymEngine::Basic;
    using SymEngine::RCP;
    using SymEngine::Symbol;

    // 定义矩形结构体用于内部快速计算
    struct Rect
    {
        double xmin, xmax, ymin, ymax;
    };

    /**
     * 将曲线点集根据矩形区域分割
     * * 输入:
     * x0, y0: 曲线坐标点
     * rects_sym: 矩形区域定义，来自 AllRegions::regions_area
     * 格式为 vector<vector<Symbol>>, 每行 [xmin, xmax, ymin, ymax]
     * * 输出 (通过引用参数):
     * x_segments, y_segments: 分割后的线段点集
     * rect_ids: 对应矩形的索引 (0-based)
     */
    inline void split_curve_by_rects_by_points(
        const std::vector<double> &x0,
        const std::vector<double> &y0,
        const std::vector<std::vector<RCP<const Symbol>>> &rects_sym,
        std::vector<std::vector<double>> &x_segments,
        std::vector<std::vector<double>> &y_segments,
        std::vector<int> &rect_ids)
    {
        // 清空输出
        x_segments.clear();
        y_segments.clear();
        rect_ids.clear();

        size_t N = x0.size();
        if (N == 0 || N != y0.size())
            return;

        // 1. 预处理：将 SymEngine 的矩形转换为 double 以提高比较效率
        // 假设 regions_area 中的 Symbol 实际上包含的是具体数值 (RealDouble)
        std::vector<Rect> rects_num;
        rects_num.reserve(rects_sym.size());

        for (const auto &row : rects_sym)
        {
            if (row.size() < 4)
                continue;
            Rect r;
            try
            {
                // eval_double 会尝试求值，如果是纯符号且没有替换值可能会报错，
                // 这里假设传入的是包含数值的 Symbol
                r.xmin = SymEngine::eval_double(*row[0]);
                r.xmax = SymEngine::eval_double(*row[1]);
                r.ymin = SymEngine::eval_double(*row[2]);
                r.ymax = SymEngine::eval_double(*row[3]);
                rects_num.push_back(r);
            }
            catch (...)
            {
                // 如果转换失败（例如纯符号），默认为无效区域或极大范围，视具体逻辑而定
                // 这里暂且设为 NaN 或跳过
                std::cerr << "[SplitCurveUtil] Warning: Failed to evaluate rectangle bounds to double." << std::endl;
                rects_num.push_back({0, 0, 0, 0});
            }
        }

        size_t seg_start = 0; // C++ 0-based index (MATLAB was 1)

        while (seg_start < N)
        {
            double cur_x = x0[seg_start];
            double cur_y = y0[seg_start];

            // 当前点属于哪些矩形
            int rect_id = -1;
            for (int i = 0; i < (int)rects_num.size(); ++i)
            {
                const auto &r = rects_num[i];
                // MATLAB: x >= xmin && x <= xmax ...
                if (cur_x >= r.xmin && cur_x <= r.xmax &&
                    cur_y >= r.ymin && cur_y <= r.ymax)
                {
                    rect_id = i;
                    break; // find(..., 1) 取第一个
                }
            }

            if (rect_id == -1)
            {
                // 如果当前点不在任何矩形中，跳到下一点
                seg_start++;
                continue;
            }

            size_t seg_end = seg_start;
            const auto &target_rect = rects_num[rect_id];

            // 找连续点仍在该矩形内
            while (seg_end + 1 < N)
            {
                double next_x = x0[seg_end + 1];
                double next_y = y0[seg_end + 1];

                if (next_x >= target_rect.xmin && next_x <= target_rect.xmax &&
                    next_y >= target_rect.ymin && next_y <= target_rect.ymax)
                {
                    seg_end++;
                }
                else
                {
                    break;
                }
            }

            // 保存该段
            std::vector<double> seg_x;
            std::vector<double> seg_y;
            seg_x.reserve(seg_end - seg_start + 1);
            seg_y.reserve(seg_end - seg_start + 1);

            for (size_t k = seg_start; k <= seg_end; ++k)
            {
                seg_x.push_back(x0[k]);
                seg_y.push_back(y0[k]);
            }

            x_segments.push_back(seg_x);
            y_segments.push_back(seg_y);
            rect_ids.push_back(rect_id); // 这里存的是 0-based 索引

            // 下一个段的起点
            seg_start = seg_end + 1;
        }
    }

} // namespace SplitCurveUtil