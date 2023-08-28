#include <chrono>
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>
#include <cstring>
#include <assert.h>
#include <unistd.h>
#include <iostream>
#include "logger.hpp"
#include "sketch.h"
#include "hash.hpp"
#include "util.h"
using namespace std;

vector<double> CountMinSketch::ground_truth(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream)
{
    seed_t* seeds = new seed_t[nrows];
    int** nt = new int*[nrows];
    for (int i=0;i<nrows;i++)
    {
        seeds[i]=clock();
        sleep(1);
        nt[i]=new int[ncols];
        memset(nt[i], 0, sizeof(int)*ncols);
    }
    seed_t sseed=clock();

    int n=stream.TOTAL_PACKETS;
    srand(clock());
    for (int i=0;i<n;i++)
    {
        data_t cur=stream.raw_data[i];
        for (int k=0;k<nrows;k++)
        {
            int pos = HASH::hash(cur, seeds[k]) % ncols;
            nt[k][pos] += 1;
        }
    }

    vector<double> rst(constraints.size(), 0);
    for (const auto& it : stream.counter)
    {
        int tpcnt=INT32_MAX;
        for (int k=0;k<nrows;k++)
        {
            int pos = HASH::hash(it.first, seeds[k]) % ncols;
            tpcnt=min(tpcnt, nt[k][pos]);
        }

        for (int i=0;i<constraints.size();i++)
        {
            if (tpcnt-it.second > constraints[i].err)
                rst[i]++;
        }
    }

    for (int i=0;i<constraints.size();i++)
        rst[i]=rst[i]/stream.getTotalFlows();
    return rst;
}

vector<double> CountMinSketch::trivial_simulate(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, int max_iter)
{
    int* nt=new int[ncols];
    seed_t seeds;
    vector<double> rst(constraints.size(), 0);
    default_random_engine gen;

    for (int u=0; u<max_iter; u++)
    {
        memset(nt, 0, ncols*sizeof(int));
        seeds=clock();
        bool stop_flag=true;

        count_t nflows=stream.getTotalFlows();
        for (count_t i=0;i<nflows;i++)
        {
            int freq=stream.new_stream();
            int pos=HASH::hash(i, seeds) % ncols;
            nt[pos]+=freq;
        }
        sort(nt, nt+ncols);

        for (int i=0;i<constraints.size();i++)
        {
            double curp = 1 - double(upper_bound(nt, nt+ncols, constraints[i].err)-nt)/ncols;
            double tpcur=(curp+rst[i]*u)/(u+1);
            if (fabs(tpcur-rst[i]) > STOP_THRESHOLD*rst[i])
            {
                stop_flag=false;
            }
            rst[i]=tpcur;
        }

        if (stop_flag)
            break;
    }
    
    for (int i=0;i<constraints.size();i++)
    {
        rst[i]=pow(rst[i], nrows);
    }
    return rst;
}

std::vector<double> CountMinSketch::simulate_without_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) 
{   
    if (!trunc)
    {
        vector<double> rst(constraints.size(), 0);
        default_random_engine gen;
        poisson_distribution<int> dist(double(stream.getTotalFlows())/ncols);
        vector<bool> stop_arr(constraints.size(), false);
        for (int u=0;;u++)
        {
            bool stop_flag=true;
            vector<double> curst(constraints.size(), 0);
            for (int i=0;i<EPOCH;i++)
            {
                count_t cur=stream.new_stream();
                count_t counter = cur;
                
                int num=dist(gen);
                for (int k=0;k<num;k++)
                {
                    count_t tp=stream.new_stream();
                    counter += tp;
                }

                int err=counter-cur;
                err=abs(err);
                for (int k=0;k<constraints.size();k++)
                {
                    if (err>constraints[k].err)
                        curst[k]++;
                }
            }

            for (int k=0;k<constraints.size();k++)
            {
                curst[k]=curst[k]/EPOCH;
                curst[k]=(curst[k]+rst[k]*u)/(u+1);

                if (!stop(u+1, rst[k], curst[k], constraints[k].prob))
                    stop_flag=false;

                rst[k]=curst[k];
            }

            if (stop_flag)
                break;
        }

        for (int i=0;i<rst.size();i++)
        {
            rst[i]=pow(rst[i], nrows);
        }
        return rst;
    }
    else
    {
        LOG_ERROR("TODO");
        exit(-1);
    }
}

std::vector<double> CountMinSketch::simulate_with_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc)
{
    double l=double(stream.getTotalFlows())/ncols;
    arcinit(ncols,constraints,stream,trunc);
    vector<double> arcprob;
    for (int i=0;i<constraints.size();i++)
    {
        int done=poisson_trunc(ncols, constraints[i], stream);
        double ppp=exp(-l);
        double s=0;
        for (int k=1;k<=done;k++)
        {
            ppp=ppp*l/k;
            s+=arctables[i][k-1]*ppp;
        }
        arcprob.push_back(pow(s,nrows));
    }
    return arcprob;
}

void CountMinSketch::arcinit(double l, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc)
{
    for (int i=0;i<arctables.size();i++)
    {
        int done=poisson_trunc(l, constraints[i], stream);
        if (arctables[i].size() >= done)
            continue;
        
        int let=arctables[i].size()+1;
        if (!trunc)
        {
            for (int n=let; n<=done; n++)
            {
                double arcnt=0;
                for (int iii=0;;iii++)
                {
                    double curcnt=0;

                    for (int v=0;v<EPOCH;v++)
                    {
                        int err=0;
                        for (int j=0;j<n;j++)
                            err+=stream.new_stream();
                        if (err>constraints[i].err)
                            curcnt++;
                    }

                    curcnt=curcnt/EPOCH;
                    curcnt=(curcnt+arcnt*iii)/(iii+1);
                    if (stop(iii+1, arcnt, curcnt))
                    {
                        arcnt=curcnt;
                        break;
                    }
                    else
                    {
                        arcnt=curcnt;
                    }
                }
                arctables[i].push_back(arcnt);
            }
        }
        else
        {
            double arcnt=0;
            truncinit(done, stream);
            default_random_engine gen;
            for (int n=let; n<=done; n++)
            {
                double arcnt=0;
                for (int iii=0;;iii++)
                {
                    double curcnt=0;
                    for (int v=0;v<EPOCH;v++)
                    {
                        int err=0;
                        int curnum=binomial_distribution<int>(n, BIG_PERCENT)(gen);

                        for (int j=0;j<curnum;j++)
                            err+=stream.trunc_not_small_stream();
                        
                        curcnt += P_bad(n-curnum, constraints[i].err-err);
                    }

                    curcnt=curcnt/EPOCH;
                    curcnt=(curcnt+arcnt*iii)/(iii+1);

                    if (stop(iii+1, arcnt, curcnt))
                    {
                        arcnt=curcnt;
                        break;
                    }
                    else
                    {
                        arcnt=curcnt;
                    }
                }
                arctables[i].push_back(arcnt);
            }
        }
    }
}

void CountMinSketch::truncinit(int done, StreamGen& stream)
{
    if (trunctables.size()>done)
        return;
    
    int start=trunctables.size();
    for (int n=start; n<=done; n++)
    {
        int* tpnt=new int[TINY_SIM_NUM];
        for (int i=0;i<TINY_SIM_NUM;i++)
        {
            tpnt[i]=0;
            for (int j=0;j<n;j++)
                tpnt[i]+=stream.trunc_tiny_stream();
        }
        sort(tpnt, tpnt+TINY_SIM_NUM);
        trunctables.push_back(tpnt);
    }
}

/**
 * @brief probability of the bad event where sum of [ntiny] tiny streams exceeds [sum] (i.e. P(S>sum))
 */
double CountMinSketch::P_bad(int ntiny, double sum)
{
    if (ntiny==0)
    {
        if (sum>0)
            return 0;
        else
            return 1;
    }
    if (sum <= 0)
        return 1;
    
    auto it=upper_bound(trunctables[ntiny], trunctables[ntiny]+TINY_SIM_NUM, sum);
    return 1 - (double(it-trunctables[ntiny])/TINY_SIM_NUM);
}

std::pair<double, double> CountMinSketch::strawman_simulate(StreamGen& stream, const std::deque<Constraint>& constraints, count_t sum) {
    std::vector<int> cols;
    std::vector<int> rows;
    for (auto i = 0; i < constraints.size(); i++) {
        double yibu = (1.0 * constraints[i].err) / sum;
        LOG_DEBUG("yibu=%lf,sum=%d",yibu,sum);
        double xige = constraints[i].prob;
        int col=exp(1)/yibu+1;
        int row=log(1.0/xige)+1;
        LOG_DEBUG("col=%d, row=%d",col,row);
        cols.push_back(col);
        rows.push_back(row);
    }
    int maxcol=*max_element(cols.begin(),cols.end());
    int maxrow=*max_element(rows.begin(),rows.end());
    return make_pair(maxcol,maxrow);
}