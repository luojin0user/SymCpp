#include "Region.hpp"
#include <iostream>

// 使用命名空间以便于访问 SymEngine 类
using SymEngine::Basic;
using SymEngine::RCP;
using SymEngine::Symbol;

Region::Region(int idx, CaseType type,
               const std::vector<RCP<const Symbol>> &area,
               int bc_type,
               const std::vector<int> &top,
               const std::vector<int> &bottom,
               const std::vector<int> &left,
               const std::vector<int> &right,
               const std::vector<int> &ES_regions_idx,
               int H_max, int N_max,
               const RCP<const Symbol> &mu_r,
               const RCP<const Basic> &J_r,
               int all_regions_num)
{
    // L/R/T/B 为 bool 值，表示边界是否存在，如果为True，说明这个边界不存在(即有邻居)
    // obj.Ln = isempty(left);
    this->Ln = left.empty();
    this->Rn = right.empty();
    this->Tn = top.empty();
    this->Bn = bottom.empty();

    this->case_type = type;
    this->idx = idx;

    // xl = area(1); ...
    // C++ vector index starts at 0
    if (area.size() < 4)
        throw std::runtime_error("Area vector must have 4 elements [xl, xr, yl, yt]");
    RCP<const Symbol> xl = area[0];
    RCP<const Symbol> xr = area[1];
    RCP<const Symbol> yl = area[2];
    RCP<const Symbol> yt = area[3];

    // Factory Pattern Implementation
    switch (type)
    {
    case CaseType::AlleyAir:
        this->impl = std::make_shared<AlleyAir>(idx, xl, xr, yl, yt, this->Ln, this->Rn, this->Tn, this->Bn, H_max, N_max, mu_r);
        break;
    case CaseType::BTAir: // NormalAir logic in MATLAB comments?
        this->impl = std::make_shared<BTAir>(idx, xl, xr, yl, yt, this->Ln, this->Rn, this->Tn, this->Bn, H_max, N_max, mu_r);
        break;
    case CaseType::NormalAir:
        this->impl = std::make_shared<NormalAir>(idx, xl, xr, yl, yt, this->Ln, this->Rn, this->Tn, this->Bn, H_max, N_max, mu_r);
        break;
    case CaseType::FerriteCurrent:
        this->impl = std::make_shared<FerriteCurrent>(idx, xl, xr, yl, yt, this->Ln, this->Rn, this->Tn, this->Bn, H_max, N_max, mu_r, J_r);
        break;
    case CaseType::Aluminum:
        // Assuming Aluminum class signature matches others
        this->impl = std::make_shared<Aluminum>(idx, xl, xr, yl, yt, this->Ln, this->Rn, this->Tn, this->Bn, H_max, N_max, mu_r);
        break;
    default:
        throw std::runtime_error("Unsupported case type");
    }

    // obj.impl.num_idx_hn = zeros(all_regions_num, 2, 'uint32');
    // 假设 BasicCase 中定义了该成员。由于它是动态添加的属性，我们在C++基类中可能没有定义。
    // 如果 BasicCase.hpp 中没有定义这些成员，这里会报错。
    // 假设在 BasicCase 中补充定义了:
    // std::vector<std::pair<int, int>> num_idx_hn; (size initialized here)
    // this->impl->num_idx_hn.resize(all_regions_num, {0, 0});

    // obj.impl.BCfuncs_loc_map = zeros(2, all_regions_num); -> std::map in our previous impl
    // this->impl->BCfuncs_loc_map.clear(); // Map naturally handles "sparse" zeros

    // obj.impl.ES_regions = false(1, all_regions_num);
    this->impl->ES_regions.assign(all_regions_num, false);

    // obj.impl.ES_regions(ES_regions) = true;
    for (int es_idx : ES_regions_idx)
    {
        if (es_idx >= 0 && es_idx < all_regions_num)
        {
            this->impl->ES_regions[es_idx] = true;
        }
    }

    this->all_regions_num = all_regions_num;

    // obj.boundarys = Boundarys(obj, all_regions_num);
    // Passing 'this' to Boundarys constructor usually requires enable_shared_from_this if we want to pass shared_ptr,
    // or just pass raw pointer/reference. Assuming Boundarys takes a reference or pointer to Region.
    this->boundarys = std::make_shared<Boundarys>(this, all_regions_num);

    this->set_boundary(top, bottom, left, right, bc_type);
}

void Region::set_boundary(const std::vector<int> &top,
                          const std::vector<int> &bottom,
                          const std::vector<int> &left,
                          const std::vector<int> &right,
                          int bc_type)
{
    // obj.boundarys.bottom = bottom; ...
    this->boundarys->bottom = bottom;
    this->boundarys->top = top;
    this->boundarys->left = left;
    this->boundarys->right = right;
    this->boundarys->bc_type = bc_type;

    // obj.impl.tops = top; ...
    // Assuming Case1/Case2/BasicCase has these members (added in previous step)
    // In C++, we need to cast if these are not in BasicCase, but we added them to Case1/2.
    // However, 'impl' is shared_ptr<BasicCase>. We should ensure BasicCase has these or cast.
    // In the previous strict C++ implementation, we put 'tops', 'bottoms' in Case1.
    // So strictly we should dynamic_cast, but assuming BasicCase defines them virtually or we added them to BasicCase header:
    // (In the previous turns, I added them to Case1.hpp. To be safe, let's assume BasicCase has them or we cast).

    // Try to downcast to Case1 (where tops/bottoms were defined in previous turn)
    if (auto case1_ptr = std::dynamic_pointer_cast<Case1>(this->impl))
    {
        case1_ptr->tops = top;
        case1_ptr->bottoms = bottom;
        case1_ptr->lefts = left;
        case1_ptr->rights = right;
    }
    // If Aluminum inherits from Case1, this works.

    // obj.all_edge_regions = false(obj.all_regions_num, 1);
    this->all_edge_regions.assign(this->all_regions_num, false);

    // regions = [top, bottom, left, right];
    std::vector<int> regions;
    regions.insert(regions.end(), top.begin(), top.end());
    regions.insert(regions.end(), bottom.begin(), bottom.end());
    regions.insert(regions.end(), left.begin(), left.end());
    regions.insert(regions.end(), right.begin(), right.end());

    // obj.all_edge_regions(regions) = true;
    for (int r_idx : regions)
    {
        if (r_idx >= 0 && r_idx < this->all_regions_num)
        {
            this->all_edge_regions[r_idx] = true;
        }
    }
}

std::vector<RCP<const Basic>> Region::get_region_solution_func()
{
    // obj.impl.gen_solution_func();
    // Using the overloaded version or default based on previous impl
    this->impl->gen_solution_func(this->Ln, this->Rn, this->Tn, this->Bn);

    // region_func = {obj.impl.eq_A_z, obj.impl.eq_B_x};
    // Need to access these fields. If they are in Case1/Case2, need cast or accessors.
    // Assuming they are promoted to BasicCase or we cast.
    // In previous Case1/2 impl, eq_A_z is public.

    std::vector<RCP<const Basic>> ret;
    if (auto c1 = std::dynamic_pointer_cast<Case1>(this->impl))
    {
        ret.push_back(c1->eq_A_z);
        ret.push_back(c1->eq_B_x);
    }
    else if (auto c2 = std::dynamic_pointer_cast<Case2>(this->impl))
    {
        ret.push_back(c2->eq_A_z);
        ret.push_back(c2->eq_B_x);
    }

    return ret;
}

Region::CoefficientResult Region::gen_region_coefficient_func(const std::vector<std::shared_ptr<Region>> &all_regions)
{
    this->all_regions = all_regions;

    // obj.impl.all_regions = all_regions;
    // We need to pass the shared_ptr<BasicCase> from the Region to the impl's all_regions list.
    // This requires converting vector<shared_ptr<Region>> to vector<shared_ptr<BasicCase>>
    // inside the impl logic, OR the impl stores weak_ptr<Region>.
    // Assuming impl has: std::vector<std::shared_ptr<BasicCase>> all_regions;
    // We update the impl's view of the world.

    std::vector<std::shared_ptr<BasicCase>> impl_regions;
    for (const auto &r : all_regions)
    {
        impl_regions.push_back(r->impl);
    }

    // Downcast to access 'all_regions' member of Case1/Case2
    if (auto c1 = std::dynamic_pointer_cast<Case1>(this->impl))
    {
        c1->all_regions = impl_regions;
    }

    // obj.boundarys.cal_BC(...)
    this->boundarys->cal_BC(this->Ln, this->Rn, this->Tn, this->Bn);

    // obj.impl.gen_coefficient_func();
    this->impl->gen_coefficient_func();

    // Collect results
    CoefficientResult result;

    // funcs = [obj.impl.eq_c0x; obj.impl.eq_c_hx; obj.impl.eq_d0x; obj.impl.eq_d_hx; obj.impl.eq_e_ny; obj.impl.eq_f_ny];
    // In our C++ impl, eq_c0x is vector<RCP>, eq_c_hx is vector<vector<RCP>> etc.
    // We need to standardize the output structure.

    // For simplicity here, assuming specific Cast to Case1 (which covers Case2 structure mostly)
    if (auto c1 = std::dynamic_pointer_cast<Case1>(this->impl))
    {
        // Since eq_c_hx is vector<vector>, and eq_c0x is vector, we push distinct groups.
        // Or we flatten everything into one huge vector.
        // MATLAB [A; B] concatenates.

        // Helper to append 2D vectors
        auto append_2d = [&](const std::vector<std::vector<RCP<const Basic>>> &src)
        {
            result.funcs.insert(result.funcs.end(), src.begin(), src.end());
        };
        // Helper to append 1D vectors (wrapped as rows)
        auto append_1d = [&](const std::vector<RCP<const Basic>> &src)
        {
            result.funcs.push_back(src);
        };

        // Case1 specific logic might differ from Case2 structure slightly in C++ class defs,
        // but here is the general idea based on MATLAB code order:

        // Note: In Case2, eq_c0x is vector<RCP>. In Case1 it might be cell(1,6) -> vector.
        // We need to handle the types carefully.

        if (auto c2 = std::dynamic_pointer_cast<Case2>(this->impl))
        {
            // Case 2 logic
            append_1d(c2->eq_c0x);
            append_1d(c2->eq_c_hx); // In Case2 impl, eq_c_hx was 1D vector
            append_1d(c2->eq_d0x);
            append_1d(c2->eq_d_hx);
            append_1d(c2->eq_e_ny);
            append_1d(c2->eq_f_ny);
        }
        else
        {
            // Case 1 logic
            // In Case1 impl, eq_c0x was not used/defined as member in provided C++ snippet,
            // but eq_c_hx was vector<vector>.
            // Using generic access or dynamic check:
            append_2d(c1->eq_c_hx);
            append_2d(c1->eq_d_hx);
            append_1d(c1->eq_e_ny);
            append_1d(c1->eq_f_ny);
        }

        result.BCfuncs_loc_map = c1->BCfuncs_loc_map;
        result.ESfuncs = c1->eq_ES;
    }

    return result;
}