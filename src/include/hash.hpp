#pragma once
#include <stdint.h>
#include "defs.hpp"

namespace HASH
{
    uint64_t hash(const data_t& data, seed_t seed = 0U);
}