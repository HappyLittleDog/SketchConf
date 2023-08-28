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

#include "include/sketch.h"

static int C(int n, int k)
{
    uint64_t ans=1;
    for (int i=0;i<k;i++)
        ans*=(n-i);
    for (int i=1;i<=k;i++)
        ans/=i;
    return ans;
}

inline int trivial_hash()
{
    if (rand()%2)
        return 1;
    else
        return -1;
}

inline int trivial_hash(data_t item, seed_t seed)
{
    if (HASH::hash(item, seed) % 2)
        return 1;
    else
        return -1;
}

vector<double> CountSketch::ground_truth(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream)
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
            nt[k][pos] += trivial_hash(cur, sseed);
        }
    }

    vector<double> rst(constraints.size(), 0);
    for (const auto& it : stream.counter)
    {
        vector<int> tpcnt;
        for (int k=0;k<nrows;k++)
        {
            int pos = HASH::hash(it.first, seeds[k]) % ncols;
            tpcnt.push_back(nt[k][pos]);
        }
        sort(tpcnt.begin(), tpcnt.end());
        
        int curst;
        if (nrows % 2)
        {
            curst=tpcnt[(nrows-1)/2];
        }
        else
        {
            curst=(tpcnt[nrows/2 - 1]+tpcnt[nrows/2]) / 2;
        }
        curst *=  trivial_hash(it.first, sseed);

        for (int i=0;i<constraints.size();i++)
        {
            if (abs(int(curst-it.second)) > constraints[i].err)
                rst[i]++;
        }
    }

    for (int i=0;i<constraints.size();i++)
        rst[i]=rst[i]/stream.getTotalFlows();
    return rst;
}

vector<double> CountSketch::trivial_simulate(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, int max_iter)
{
    if (nrows % 2 == 0)
    {
        LOG_DEBUG("Unsupported");
        return vector<double>(constraints.size(), 1);
    }
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
            nt[pos]+=freq * trivial_hash();
        }
        sort(nt, nt+ncols);

        for (int i=0;i<constraints.size();i++)
        {
            double curp=1 - double(upper_bound(nt, nt+ncols, constraints[i].err)-nt)/ncols;
            int st=(nrows+1)/2;
            double tpans=0;
            for (int k=st; k<=nrows; k++)
            {
                tpans+=C(nrows, k)*pow(curp, k)*pow(1-curp, nrows-k);
            }
            tpans *= 2;

            double tpcur=(tpans+rst[i]*u)/(u+1);
            if (fabs(tpcur-rst[i]) > STOP_THRESHOLD*rst[i])
            {
                stop_flag=false;
            }
            rst[i]=tpcur;
        }

        if (stop_flag)
            break;
    }
    
    return rst;
}

vector<double> CountSketch::simulate_without_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc) 
{
    if (nrows % 2 == 0)
    {
        LOG_DEBUG("Unsupported");
        return vector<double>(constraints.size(), 1);
    }
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
                int counter = 0;
                
                int num=dist(gen);
                for (int k=0;k<num;k++)
                {
                    int tp=stream.new_stream() * trivial_hash();
                    counter += tp;
                }
                
                for (int k=0;k<constraints.size();k++)
                {
                    if (counter>constraints[k].err)
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
            {
                break;
            }
        }

        int st=(nrows+1)/2;
        for (int i=0;i<rst.size();i++)
        {
            double ans=0;
            for (int k=st; k<=nrows; k++)
            {
                ans+=C(nrows, k)*pow(rst[i], k)*pow(1-rst[i], nrows-k);
            }
            rst[i]=ans*2;
        }
        return rst;
    }
    else 
    {
        LOG_ERROR("TODO");
        exit(-1);
    }
}

vector<double> CountSketch::simulate_with_reuse(count_t nrows, count_t ncols, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc)
{
    if (nrows % 2 == 0)
    {
        LOG_DEBUG("Unsupported");
        return vector<double>(constraints.size(), 1);
    }

    double l;
    double done;
    int nflows=stream.getTotalFlows();
    l = double(nflows)/ncols;
    arcinit(ncols,constraints,stream,trunc);
    vector<double> arcprob;
    int st=(nrows+1)/2;
    for (int u=0;u<constraints.size();u++)
    {
        done=poisson_trunc(ncols, constraints[u], stream);

        double ppp=exp(-l);
        double s=ppp*arctables[u][0];
        for (int k=1;k<=done;k++)
        {
            ppp=ppp*l/k;
            s+=arctables[u][k]*ppp;
        }

        double ans=0;
        for (int k=st; k<=nrows; k++)
        {
            ans+=C(nrows, k)*pow(s, k)*pow(1-s, nrows-k);
        }
        arcprob.push_back(ans*2);
    }
    return arcprob;
}

void CountSketch::arcinit(double l, const std::deque<Constraint>& constraints, StreamGen& stream, bool trunc)
{
    for (int i=0;i<arctables.size();i++)
    {
        int done=poisson_trunc(l, constraints[i], stream);
        if (arctables[i].size() > done)
            continue;
        
        int let=arctables[i].size();
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
                            err+=stream.new_stream()*trivial_hash();
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
                            err+=stream.trunc_not_small_stream() * trivial_hash();
                        
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

void CountSketch::truncinit(int done, StreamGen& stream)
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
                tpnt[i]+=stream.trunc_tiny_stream() * trivial_hash();
        }
        sort(tpnt, tpnt+TINY_SIM_NUM);
        trunctables.push_back(tpnt);
    }
}

double CountSketch::P_bad(int ntiny, double sum)
{
    if (ntiny==0)
    {
        if (sum>0)
            return 0;
        else
            return 1;
    }
    
    auto it=upper_bound(trunctables[ntiny], trunctables[ntiny]+TINY_SIM_NUM, sum);
    return 1 - (double(it-trunctables[ntiny])/TINY_SIM_NUM);
}

pair<double,double> CountSketch::strawman_simulate(StreamGen& stream, const std::deque<Constraint>& constraints,count_t sum)
{
    vector<double> ds;
    vector<double> ws;
    double l2=stream.L2_Norm();
    for (auto& it : constraints)
    {
        ds.push_back(ceil(8*log(1/it.prob)));
        ws.push_back(ceil(4*l2*l2/it.err/it.err));
    }
    sort(ds.begin(), ds.end());
    sort(ws.begin(), ws.end());
    return make_pair(ws.back(), ds.back());
}