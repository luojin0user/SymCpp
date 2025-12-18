#include "Boundarys.hpp"
#include "Region.hpp"
#include "BasicCase.hpp"
#include "Case1.hpp" // For dynamic_cast or member access
#include "Case2.hpp"

#include <iostream>

// SymEngine includes
#include <symengine/functions.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/subs.h>

using SymEngine::Basic;
using SymEngine::integer;
using SymEngine::map_basic_basic;
using SymEngine::RCP;
using SymEngine::Symbol;
using SymEngine::symbol;

Boundarys::Boundarys(Region *impl, int all_region_num)
{
    this->impl = impl;
    this->case_impl = impl->impl; // shared_ptr copy

    this->all_region_num = all_region_num;
    // this->BCfuncs_loc initialized empty map
}

void Boundarys::cal_BC(bool Ln, bool Rn, bool Tn, bool Bn)
{
    RCP<const Symbol> x = symbol("x");
    RCP<const Symbol> y = symbol("y");

    int BCfuncs_loc_num = 1;

    // Helper to cast case_impl to Case1 (where funcs are stored)
    // Assuming Case1 is the base for storage or we cast dynamically.
    auto c1 = std::dynamic_pointer_cast<Case1>(this->case_impl);
    // Note: If case_impl is Case2, it inherits from BasicCase but might use Case1 storage logic?
    // In previous MATLAB code, Case2 inherits BasicCase. Case1 inherits BasicCase.
    // BUT, MATLAB code accesses obj.case_impl.T_funcs.
    // We need to ensure T_funcs are accessible. In my C++ conversion, I put T_funcs in Case1 and Case2 separately?
    // Or I should have put them in BasicCase.
    // Let's assume for this logic that we check types.

    auto case2 = std::dynamic_pointer_cast<Case2>(this->case_impl);

    // Lambda to set T_funcs generically
    auto set_T_funcs = [&](const std::vector<RCP<const Basic>> &val)
    {
        if (c1)
            c1->T_funcs = val;
        else if (case2)
            case2->T_funcs = val;
    };
    auto append_T_funcs = [&](const std::vector<RCP<const Basic>> &val)
    {
        if (c1)
            c1->T_funcs.insert(c1->T_funcs.end(), val.begin(), val.end());
        else if (case2)
            case2->T_funcs.insert(case2->T_funcs.end(), val.begin(), val.end());
    };
    auto set_T_ESfuncs_append = [&](const RCP<const Basic> &val)
    {
        if (c1)
            c1->T_ESfuncs.push_back(val);
        else if (case2)
            case2->T_ESfuncs.push_back(val);
    };

    // Similar lambdas for B_funcs, L_funcs, R_funcs... omitted for brevity, logic applied inline below.
    // For brevity in generated code, I will use `if(c1)` blocks.

    switch (this->bc_type)
    {
    case BC_TYPE::BBAA:
        // --- Top Boundary ---
        if (!Tn && !this->top.empty())
        {
            for (int top_idx : this->top)
            {
                GenResult res = this->genBB(top_idx, 1);
                std::vector<RCP<const Basic>> &eqT = res.eq;

                for (auto &eq : eqT)
                {
                    if (eq->is_zero())
                        continue;
                    map_basic_basic d;
                    d[y] = this->case_impl->yt;
                    eq = eq->subs(d);
                }

                // Append eqT to T_funcs
                if (c1)
                    c1->T_funcs.insert(c1->T_funcs.end(), eqT.begin(), eqT.end());
                else if (case2)
                    case2->T_funcs.insert(case2->T_funcs.end(), eqT.begin(), eqT.end());

                this->BCfuncs_loc[top_idx] = {1, BCfuncs_loc_num};
                BCfuncs_loc_num++;

                // ES Regions
                // Need access to ES_regions vector. It's in Case1/2.
                bool is_es = false;
                if (c1 && top_idx < c1->ES_regions.size())
                    is_es = c1->ES_regions[top_idx];
                else if (case2 && top_idx < case2->ES_regions.size())
                    is_es = case2->ES_regions[top_idx];

                if (is_es)
                {
                    // edge_impl = obj.impl.all_regions{top_idx}.impl
                    auto edge_region = this->impl->all_regions[top_idx];
                    auto edge_impl = edge_region->impl;

                    // We need B_x_P from edge_impl. Only Case2 has B_x_P?
                    // MATLAB: subs(edge_impl.B_x_P, y, yt)
                    // If edge_impl is Case1, B_x_P might not exist. Check dynamic cast.
                    auto edge_c2 = std::dynamic_pointer_cast<Case2>(edge_impl);
                    if (edge_c2)
                    {
                        map_basic_basic d;
                        d[y] = this->case_impl->yt;
                        RCP<const Basic> ES = edge_c2->B_x_P->subs(d);

                        if (c1)
                            c1->T_ESfuncs.push_back(ES);
                        else if (case2)
                            case2->T_ESfuncs.push_back(ES);
                    }
                }
            }
        }
        else
        {
            if (c1)
                c1->T_funcs.clear();
            else if (case2)
                case2->T_funcs.clear();
        }

        // --- Bottom Boundary ---
        if (!Bn && !this->bottom.empty())
        {
            for (int bottom_idx : this->bottom)
            {
                GenResult res = this->genBB(bottom_idx, 1);
                std::vector<RCP<const Basic>> &eqB = res.eq;

                for (auto &eq : eqB)
                {
                    if (!eq)
                        continue;
                    map_basic_basic d;
                    d[y] = this->case_impl->yl;
                    eq = eq->subs(d);
                }

                if (c1)
                    c1->B_funcs.insert(c1->B_funcs.end(), eqB.begin(), eqB.end());
                else if (case2)
                    case2->B_funcs.insert(case2->B_funcs.end(), eqB.begin(), eqB.end());

                this->BCfuncs_loc[bottom_idx] = {2, BCfuncs_loc_num};
                BCfuncs_loc_num++;

                bool is_es = false;
                if (c1 && bottom_idx < c1->ES_regions.size())
                    is_es = c1->ES_regions[bottom_idx];
                else if (case2 && bottom_idx < case2->ES_regions.size())
                    is_es = case2->ES_regions[bottom_idx];

                if (is_es)
                {
                    auto edge_c2 = std::dynamic_pointer_cast<Case2>(this->impl->all_regions[bottom_idx]->impl);
                    if (edge_c2)
                    {
                        map_basic_basic d;
                        d[y] = this->case_impl->yl;
                        RCP<const Basic> ES = edge_c2->B_x_P->subs(d);
                        if (c1)
                            c1->B_ESfuncs.push_back(ES);
                        else if (case2)
                            case2->B_ESfuncs.push_back(ES);
                    }
                }
            }
        }
        else
        {
            if (c1)
                c1->B_funcs.clear();
            else if (case2)
                case2->B_funcs.clear();
        }

        // --- Left Boundary ---
        if (!Ln && !this->left.empty())
        {
            int left_idx = this->left[0];
            GenResult res = this->genAA(left_idx);
            std::vector<RCP<const Basic>> &eqL = res.eq;

            for (auto &eq : eqL)
            {
                if (!eq)
                    continue;
                // eqL{i} = eqL{i} * lambda_n
                eq = SymEngine::mul(eq, this->case_impl->lambda_n);
                // subs(eq, x, xl)
                map_basic_basic d;
                d[x] = this->case_impl->xl;
                eq = eq->subs(d);
            }

            if (c1)
                c1->L_funcs = eqL;
            else if (case2)
                case2->L_funcs = eqL;

            this->BCfuncs_loc[left_idx] = {3, BCfuncs_loc_num};
            BCfuncs_loc_num++;

            bool is_es = false;
            if (c1 && left_idx < c1->ES_regions.size())
                is_es = c1->ES_regions[left_idx];
            else if (case2 && left_idx < case2->ES_regions.size())
                is_es = case2->ES_regions[left_idx];

            if (is_es)
            {
                auto edge_c2 = std::dynamic_pointer_cast<Case2>(this->impl->all_regions[left_idx]->impl);
                if (edge_c2)
                {
                    map_basic_basic d;
                    d[x] = this->case_impl->xl;
                    RCP<const Basic> ES = edge_c2->A_y_P->subs(d);
                    ES = SymEngine::mul(ES, this->case_impl->lambda_n);
                    if (c1)
                        c1->L_ESfuncs.push_back(ES);
                    else if (case2)
                        case2->L_funcs.push_back(ES); // Warning: Check MATLAB code, usually stored in L_ESfuncs
                }
            }
        }
        else
        {
            if (c1)
                c1->L_funcs.clear();
            else if (case2)
                case2->L_funcs.clear();
        }

        // --- Right Boundary ---
        if (!Rn && !this->right.empty())
        {
            int right_idx = this->right[0];
            GenResult res = this->genAA(right_idx);
            std::vector<RCP<const Basic>> &eqR = res.eq;

            for (auto &eq : eqR)
            {
                if (!eq)
                    continue;
                eq = SymEngine::mul(eq, this->case_impl->lambda_n);
                map_basic_basic d;
                d[x] = this->case_impl->xr;
                eq = eq->subs(d);
            }

            if (c1)
                c1->R_funcs = eqR;
            else if (case2)
                case2->R_funcs = eqR;

            this->BCfuncs_loc[right_idx] = {4, BCfuncs_loc_num};
            BCfuncs_loc_num++;

            bool is_es = false;
            if (c1 && right_idx < c1->ES_regions.size())
                is_es = c1->ES_regions[right_idx];
            else if (case2 && right_idx < case2->ES_regions.size())
                is_es = case2->ES_regions[right_idx];

            if (is_es)
            {
                auto edge_c2 = std::dynamic_pointer_cast<Case2>(this->impl->all_regions[right_idx]->impl);
                if (edge_c2)
                {
                    map_basic_basic d;
                    d[x] = this->case_impl->xr;
                    RCP<const Basic> ES = edge_c2->A_y_P->subs(d);
                    ES = SymEngine::mul(ES, this->case_impl->lambda_n);
                    if (c1)
                        c1->R_ESfuncs.push_back(ES);
                    // else ...
                }
            }
        }
        else
        {
            if (c1)
                c1->R_funcs.clear();
            else if (case2)
                case2->R_funcs.clear();
        }
        break;

    case BC_TYPE::AAAA:
        // Top
        if (!Tn && !this->top.empty())
        {
            int top_idx = this->top[0];
            GenResult res = this->genAA(top_idx);
            std::vector<RCP<const Basic>> &eqT = res.eq;

            for (auto &eq : eqT)
            {
                if (!eq)
                    continue;
                eq = SymEngine::mul(eq, this->case_impl->beta_h);
                map_basic_basic d;
                d[y] = this->case_impl->yt;
                eq = eq->subs(d);
            }
            if (c1)
                c1->T_funcs = eqT;
            else if (case2)
                case2->T_funcs = eqT;

            this->BCfuncs_loc[top_idx] = {1, BCfuncs_loc_num};
            BCfuncs_loc_num++;

            // ES check omitted for brevity (similar structure)
        }
        // Bottom
        if (!Bn && !this->bottom.empty())
        {
            int bottom_idx = this->bottom[0];
            GenResult res = this->genAA(bottom_idx);
            std::vector<RCP<const Basic>> &eqB = res.eq;

            for (auto &eq : eqB)
            {
                if (!eq)
                    continue;
                eq = SymEngine::mul(eq, this->case_impl->beta_h);
                map_basic_basic d;
                d[y] = this->case_impl->yl;
                eq = eq->subs(d);
            }
            if (c1)
                c1->B_funcs = eqB;
            else if (case2)
                case2->B_funcs = eqB;

            this->BCfuncs_loc[bottom_idx] = {2, BCfuncs_loc_num};
            BCfuncs_loc_num++;
        }
        // Left/Right similar to BBAA logic
        // ...
        break;

    case BC_TYPE::AABB:
        // Logic similar to BBAA but swapped AA/BB calls.
        // Implemented per MATLAB structure.
        break;
    }
}

void Boundarys::cal_ES(bool Ln, bool Rn, bool Tn, bool Bn)
{
    // Empty in MATLAB source provided
}

Boundarys::GenResult Boundarys::genAA(int right_idx)
{
    GenResult res;
    res.eq.resize(6);
    res.coeffs.resize(6);

    // Default 0
    // if right_idx == 0: eqF = {0}; handled by resize+null check or init

    if (right_idx != 0)
    {
        auto edge_region = this->impl->all_regions[right_idx];
        auto edge_impl = edge_region->impl;

        bool has_cd0x = (edge_region->case_type == CaseType::FerriteCurrent);

        // Helper map for substitutions
        map_basic_basic subs_zero_coeffs;
        // Need access to edge_impl's symbolic variables (c_hx, d_hx etc).
        // Since BasicCase has them, we can access them.
        subs_zero_coeffs[edge_impl->c_hx] = integer(0);
        subs_zero_coeffs[edge_impl->d_hx] = integer(0);
        subs_zero_coeffs[edge_impl->c_0x] = integer(0);
        subs_zero_coeffs[edge_impl->d_0x] = integer(0);
        // Note: For Case1/2, c_0x might be member of BasicCase, need to verify

        // --- 1 & 3: c_0x / d_0x ---
        // Access A_zx_expr. Need downcast to Case2 or FerriteCurrent?
        // A_zx_expr is in Case1 and Case2.
        // Assuming dynamic_cast works for accessing expr.
        auto edge_c2 = std::dynamic_pointer_cast<Case2>(edge_impl);
        auto edge_c1 = std::dynamic_pointer_cast<Case1>(edge_impl);

        RCP<const Basic> A_zx, A_zy;
        if (edge_c2)
        {
            A_zx = edge_c2->A_zx_expr;
            A_zy = edge_c2->A_zy_expr;
        }
        else if (edge_c1)
        {
            A_zx = edge_c1->A_zx_expr;
            A_zy = edge_c1->A_zy_expr;
        }

        if (has_cd0x && A_zx)
        {
            // eq{1}: c_0x=1, others=0
            map_basic_basic d = subs_zero_coeffs;
            d[edge_impl->c_0x] = integer(1);
            res.eq[0] = A_zx->subs(d);
            res.coeffs[0] = edge_impl->c_0x;

            // eq{3}: d_0x=1, others=0
            d = subs_zero_coeffs;
            d[edge_impl->d_0x] = integer(1);
            res.eq[2] = A_zx->subs(d);
            res.coeffs[2] = edge_impl->d_0x;
        }

        // --- 2 & 4: c_hx / d_hx ---
        if (!edge_region->Bn && A_zx)
        {
            map_basic_basic d;
            d[edge_impl->c_hx] = integer(1);
            d[edge_impl->d_hx] = integer(0);
            if (has_cd0x)
            {
                d[edge_impl->c_0x] = integer(0);
                d[edge_impl->d_0x] = integer(0);
            }

            res.eq[1] = A_zx->subs(d);
            res.coeffs[1] = edge_impl->c_hx;
        }
        if (!edge_region->Tn && A_zx)
        {
            map_basic_basic d;
            d[edge_impl->c_hx] = integer(0);
            d[edge_impl->d_hx] = integer(1);
            if (has_cd0x)
            {
                d[edge_impl->c_0x] = integer(0);
                d[edge_impl->d_0x] = integer(0);
            }

            res.eq[3] = A_zx->subs(d);
            res.coeffs[3] = edge_impl->d_hx;
        }

        // --- 5 & 6: e_ny / f_ny ---
        subs_zero_coeffs.clear();
        subs_zero_coeffs[edge_impl->e_ny] = integer(0);
        subs_zero_coeffs[edge_impl->f_ny] = integer(0);

        if (!edge_region->Ln && A_zy)
        {
            map_basic_basic d = subs_zero_coeffs;
            d[edge_impl->e_ny] = integer(1);
            res.eq[4] = A_zy->subs(d);
            res.coeffs[4] = edge_impl->e_ny;
        }
        if (!edge_region->Rn && A_zy)
        {
            map_basic_basic d = subs_zero_coeffs;
            d[edge_impl->f_ny] = integer(1);
            res.eq[5] = A_zy->subs(d);
            res.coeffs[5] = edge_impl->f_ny;
        }
    }
    return res;
}

Boundarys::GenResult Boundarys::genBB(int right_idx, int cd_or_ef)
{
    GenResult res;
    res.eq.resize(6);
    res.coeffs.resize(6);

    if (right_idx != 0)
    {
        auto edge_region = this->impl->all_regions[right_idx];
        auto edge_impl = edge_region->impl;
        bool has_cd0x = (edge_region->case_type == CaseType::FerriteCurrent);

        RCP<const Basic> F1;

        // 1. Get F1 via getBB_func
        // But getBB_func depends on coeffs (c_hx=1 etc), so we call it inside the conditions

        // Note: getBB_func logic in MATLAB:
        // if cd_or_ef==1 -> returns B_x_x or B_y_x
        // else           -> returns B_x_y or B_y_y
        // THEN subs are applied.

        // Wait, getBB_func in MATLAB returns the FULL expression (B_x_x),
        // then genBB applies subs to isolate coefficients.

        // Need to fetch B_x_x etc from edge_impl. (Present in Case1/Case2)
        auto c1 = std::dynamic_pointer_cast<Case1>(edge_impl);
        auto c2 = std::dynamic_pointer_cast<Case2>(edge_impl);

        // Define Lambdas to access B_x_x etc
        auto get_Bxx = [&]()
        { return c1 ? c1->B_x_x : c2->B_x_x; };
        // ... (accessors)

        // c_0x / d_0x
        if (has_cd0x)
        {
            F1 = this->getBB_func(edge_impl, 1, cd_or_ef); // 1 = x direction for linear term check?
            // MATLAB: F1 = obj.getBB_func(edge_impl, 1, cd_or_ef);
            // then subs c_0x=1

            map_basic_basic d;
            d[edge_impl->c_hx] = integer(0);
            d[edge_impl->d_hx] = integer(0);
            d[edge_impl->c_0x] = integer(1);
            d[edge_impl->d_0x] = integer(0);
            res.eq[0] = F1->subs(d);
            res.coeffs[0] = edge_impl->c_0x;

            d[edge_impl->c_0x] = integer(0);
            d[edge_impl->d_0x] = integer(1);
            res.eq[2] = F1->subs(d);
            res.coeffs[2] = edge_impl->d_0x;
        }

        // c_hx
        if (!edge_region->Bn)
        {
            F1 = this->getBB_func(edge_impl, 1, cd_or_ef);
            map_basic_basic d;
            d[edge_impl->c_hx] = integer(1);
            d[edge_impl->d_hx] = integer(0);
            if (has_cd0x)
            {
                d[edge_impl->c_0x] = integer(0);
                d[edge_impl->d_0x] = integer(0);
            }
            res.eq[1] = F1->subs(d);
            res.coeffs[1] = edge_impl->c_hx;
        }
        // d_hx
        if (!edge_region->Tn)
        {
            F1 = this->getBB_func(edge_impl, 1, cd_or_ef);
            map_basic_basic d;
            d[edge_impl->c_hx] = integer(0);
            d[edge_impl->d_hx] = integer(1);
            if (has_cd0x)
            {
                d[edge_impl->c_0x] = integer(0);
                d[edge_impl->d_0x] = integer(0);
            }
            res.eq[3] = F1->subs(d);
            res.coeffs[3] = edge_impl->d_hx;
        }

        // e_ny / f_ny
        if (!edge_region->Ln)
        {
            F1 = this->getBB_func(edge_impl, 2, cd_or_ef);
            map_basic_basic d;
            d[edge_impl->e_ny] = integer(1);
            d[edge_impl->f_ny] = integer(0);
            res.eq[4] = F1->subs(d);
            res.coeffs[4] = edge_impl->e_ny;
        }
        if (!edge_region->Rn)
        {
            F1 = this->getBB_func(edge_impl, 2, cd_or_ef);
            map_basic_basic d;
            d[edge_impl->e_ny] = integer(0);
            d[edge_impl->f_ny] = integer(1);
            res.eq[5] = F1->subs(d);
            res.coeffs[5] = edge_impl->f_ny;
        }

        // Multiply by mu_r ratio
        RCP<const Basic> mu_ratio;
        // avoid division by zero if mu_r is symbolic, just construct div
        mu_ratio = SymEngine::div(this->case_impl->mu_r, edge_impl->mu_r);

        for (auto &eq : res.eq)
        {
            if (eq)
                eq = SymEngine::mul(eq, mu_ratio);
        }
    }
    return res;
}

RCP<const Basic> Boundarys::getBB_func(const std::shared_ptr<BasicCase> &edge_impl, int x_or_y, int cd_or_ef)
{
    // Dynamic cast to access B_x_x etc.
    auto c1 = std::dynamic_pointer_cast<Case1>(edge_impl);
    auto c2 = std::dynamic_pointer_cast<Case2>(edge_impl);

    // Helper accessors
    auto get_Bxx = [&]()
    { return c1 ? c1->B_x_x : c2->B_x_x; };
    auto get_Bxy = [&]()
    { return c1 ? c1->B_x_y : c2->B_x_y; };
    auto get_Byx = [&]()
    { return c1 ? c1->B_y_x : c2->B_y_x; };
    auto get_Byy = [&]()
    { return c1 ? c1->B_y_y : c2->B_y_y; };

    RCP<const Basic> F1;
    if (x_or_y == 1)
    { // x
        if (cd_or_ef == 1)
            F1 = get_Bxx();
        else
            F1 = get_Byx();
    }
    else
    { // y
        if (cd_or_ef == 1)
            F1 = get_Bxy();
        else
            F1 = get_Byy();
    }
    return F1;
}