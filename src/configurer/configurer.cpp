#include <deque>
#include <iostream>
#include <memory>
#include "include/configurer.h"
#include "include/defs.hpp"
#include "include/sketch.h"
#include "include/streamgen.h"
using namespace std;

bool satisfy(int ncols, int nrows, const Demand& demand, StreamGen& stream, bool trivial, bool trunc, bool reuse, BaseSketch* sketch) 
{
    auto errate = sketch->simulate(nrows, ncols, demand.constraints, stream, trivial, trunc, reuse);
    for (int i = 0; i < errate.size(); i++) 
    {
        if (errate[i] > demand.constraints[i].prob)
            return false;
    }
    return true;
}
pair<int, int> search(StreamGen& stream, const Demand& demand, bool trivial, bool trunc, bool reuse) 
{
    LOG_DEBUG("@search");
    std::unique_ptr<BaseSketch> sketch;
    auto type = demand.type;
    switch (type) 
    {
        case SketchType::CM:
            sketch = std::make_unique<CountMinSketch>(demand.constraints.size());
            break;
        
        case SketchType::COUNT:
            sketch = std::make_unique<CountSketch>(demand.constraints.size());
            break;
        
        case SketchType::TOWER:
            sketch = std::make_unique<TowerSketch>(demand.constraints.size());
            break;
        
        case SketchType::NITROCM:
            sketch = std::make_unique<NitroCMSketch>(demand.constraints.size(), demand.sample_rate);
            break;

        case SketchType::NITROCOUNT:
            sketch = std::make_unique<NitroCountSketch>(demand.constraints.size(), demand.sample_rate);
            break;
        
        default:
            LOG_ERROR("Unrecognized sketch type: %d", type);
    }
    int ncols = 0, nrows = 0;
    int left = 0, right = MAX_MEMORY;
    while (left + 1 < right) 
    {
        int cur = (left + right) / 8;
        cur = cur * 4;
        bool satisfied = false;
        LOG_DEBUG("search::current memory = %d", cur);
        for (int i = 1; i < 7; i++) 
        {
            LOG_DEBUG("\tcurrent col = %d, row = %d", cur / (4 * i), i);
            if (satisfy(cur / (4 * i), i, demand, stream, trivial, trunc, reuse, sketch.get())) 
            {
                LOG_DEBUG("\t\tSatisfy!");
                ncols = cur / (4 * i);
                nrows = i;
                right = cur - 4;
                satisfied = true;
                break;
            } 
            else 
            {
                LOG_DEBUG("\t\tNot satisfy...");
            }
        }
        if (!satisfied) 
        {
            left = cur + 4;
        }
    }
    LOG_RESULT("type = %d, ncols = %d, nrows = %d, total memory = %d", demand.type, ncols, nrows, ncols * nrows * 4);
    return make_pair(nrows, ncols);
}

pair<int, int> StrawmanSearch(StreamGen& stream, const Demand& demand) 
{
    auto type = demand.type;
    std::unique_ptr<BaseSketch> sketch;

    switch (type) 
    {
        case SketchType::CM:
            sketch = std::make_unique<CountMinSketch>(1);
            break;
        
        case SketchType::COUNT:
            sketch = std::make_unique<CountSketch>(1);
            break;
    
        default:
            LOG_ERROR("Unrecognized sketch type: %d", type);
    }
    double ncols, nrows;
    std::pair<double, double> ans = sketch->strawman_simulate(stream, demand.constraints, stream.TOTAL_PACKETS);
    ncols = ans.first;
    nrows = ans.second;
    LOG_RESULT("type = %d, ncols = %lf, nrows = %lf, total memory = %lf", demand.type, ncols, nrows, 4*nrows*ncols);
    return make_pair(nrows, ncols);
}
