#include "farm.hpp"
#include "hash.hpp"
#include "defs.hpp"

namespace HASH
{
    uint64_t hash(const data_t& data, seed_t seed)
    {
        return NAMESPACE_FOR_HASH_FUNCTIONS::Hash64WithSeed(reinterpret_cast<const char *>(&data), sizeof(data_t), seed);
    }
}