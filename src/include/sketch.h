#pragma once
#include "defs.hpp"
#include "streamgen.h"
#include <vector>
#include <cstring>
using namespace std;

class BaseSketch 
{
private:
    virtual std::vector<double> simulate_without_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) = 0; 
    virtual std::vector<double> simulate_with_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) = 0;  
public:
    BaseSketch() = default;
    virtual ~BaseSketch() = default;
    std::vector<double> simulate(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trivial, bool trunc, bool reuse)
    {
        if (trivial)
            return trivial_simulate(nrows, ncols, constraints, stream);
        else if (reuse)
            return simulate_with_reuse(nrows,ncols,constraints,stream,trunc);
        else
            return simulate_without_reuse(nrows,ncols,constraints,stream,trunc);
    }
    virtual std::vector<double> ground_truth(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream) {LOG_ERROR("TODO"); return std::vector<double>();};
    virtual std::pair<double,double> strawman_simulate(StreamGen& stream, const std::deque<Constraint>& constraints,count_t sum) = 0;
    virtual std::vector<double> trivial_simulate(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, int max_iter=INT32_MAX) { LOG_ERROR("TODO"); return std::vector<double>();};
};

class CountMinSketch : public BaseSketch 
{
public:
    CountMinSketch(int n_constraints)
    {
        for (int i=0; i<n_constraints; i++)
        {
            arctables.push_back(vector<double>());
        }

        int* tp=new int[TINY_SIM_NUM];
        memset(tp, 0, TINY_SIM_NUM*sizeof(int));
        trunctables.push_back(tp);
    }
private:
    vector<vector<double> > arctables;
    vector<int*> trunctables;
    virtual std::vector<double> simulate_without_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    virtual std::vector<double> simulate_with_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    void arcinit(double l, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc);
    void truncinit(int done, StreamGen& stream);
    double P_bad(int ntiny, double sum);
public:
    virtual std::vector<double> ground_truth(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream) override;
    virtual std::pair<double,double> strawman_simulate(StreamGen& stream, const std::deque<Constraint>& constraints,count_t sum) override;
    virtual std::vector<double> trivial_simulate(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, int max_iter=INT32_MAX) override;
};

class CountSketch : public BaseSketch 
{
public:
    CountSketch(int n_constraints)
    {
        for (int i=0; i<n_constraints; i++)
        {
            arctables.push_back(vector<double>());
        }

        int* tp=new int[TINY_SIM_NUM];
        memset(tp, 0, TINY_SIM_NUM*sizeof(int));
        trunctables.push_back(tp);
    }
private:
    vector<vector<double> > arctables;
    vector<int*> trunctables;
    virtual std::vector<double> simulate_without_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    virtual std::vector<double> simulate_with_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    void arcinit(double l, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc);
    void truncinit(int done, StreamGen& stream);
    double P_bad(int ntiny, double sum);
public:
    virtual std::vector<double> ground_truth(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream) override;
    virtual std::pair<double,double> strawman_simulate(StreamGen& stream, const std::deque<Constraint>& constraints,count_t sum) override;
    virtual std::vector<double> trivial_simulate(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, int max_iter=INT32_MAX) override;
};

class TowerSketch : public BaseSketch 
{
public:
    TowerSketch(int n_constraints)
    {
        for (int i=0;i<4;i++)
        {
            for (int j=0;j<n_constraints;j++)
                arctables[i].push_back(vector<double>());
        }

        int* tp=new int[TINY_SIM_NUM];
        memset(tp, 0, TINY_SIM_NUM*sizeof(int));
        trunctables.push_back(tp);
    }
private:
    vector<vector<double> > arctables[4];
    vector<int*> trunctables;
    virtual std::vector<double> simulate_without_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    virtual std::vector<double> simulate_with_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    void arcinit(vector<double >& l, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc);
    void truncinit(int done, StreamGen& stream);
    double P_bad(int ntiny, double sum);
public:
    virtual std::vector<double> ground_truth(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream) override;
    virtual std::pair<double,double> strawman_simulate(StreamGen& stream, const std::deque<Constraint>& constraints,count_t sum) override;
    virtual std::vector<double> trivial_simulate(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, int max_iter=INT32_MAX) override;
};

class NitroCMSketch : public BaseSketch 
{
public:
    NitroCMSketch(int n_constraints, double rate) : sample_rate(rate)
    {
        for (int i=0; i<n_constraints; i++)
        {
            arctables.push_back(vector<double>());
        }

        int* tp=new int[TINY_SIM_NUM];
        memset(tp, 0, TINY_SIM_NUM*sizeof(int));
        trunctables.push_back(tp);
    }
private:
    vector<vector<double> > arctables;
    vector<int*> trunctables;
    double sample_rate;
    virtual std::vector<double> simulate_without_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    virtual std::vector<double> simulate_with_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    void arcinit(double l, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc);
    void truncinit(int done, StreamGen& stream);
    double P_bad(int ntiny, double sum);
public:
    virtual std::vector<double> ground_truth(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream) override;
    virtual std::pair<double,double> strawman_simulate(StreamGen& stream, const std::deque<Constraint>& constraints,count_t sum) override;
    virtual std::vector<double> trivial_simulate(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, int max_iter=INT32_MAX) override;
};

class NitroCountSketch : public BaseSketch 
{
public:
    NitroCountSketch(int n_constraints, double rate) : sample_rate(rate)
    {
        for (int i=0; i<n_constraints; i++)
        {
            arctables.push_back(vector<double>());
        }
    }
private:
    vector<vector<double> > arctables;
    double sample_rate;
    virtual std::vector<double> simulate_without_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    virtual std::vector<double> simulate_with_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) override;
    void arcinit(double l, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc);
public:
    virtual std::vector<double> ground_truth(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream) override;
    virtual std::pair<double,double> strawman_simulate(StreamGen& stream, const std::deque<Constraint>& constraints,count_t sum) override;
    virtual std::vector<double> trivial_simulate(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, int max_iter=INT32_MAX) override;
};