#pragma once

#include <AlleyAir.hpp>
#include <Aluminum.hpp>
#include <BTAir.hpp>
#include <FerriteCurrent.hpp>
#include <NormalAir.hpp>

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
