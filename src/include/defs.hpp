#pragma once
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <deque>
#include <cstdint>
#include <string>
#include <chrono>
#include <cmath>
#include "logger.hpp"

typedef uint32_t data_t;
typedef uint32_t count_t;
typedef uint64_t seed_t;

static const double EPOCH = 1000;
static const int64_t MAX_MEMORY=20'000'000;
extern count_t TINY_THRESHOLD;
extern double BIG_PERCENT;
int TINY_SIM_NUM = 2000;
static const double STOP_THRESHOLD=1e-4;
static const double _K_=2654;
static const double PHI=2.57583; // delta = 1%
static const double RELATIVE_ERROR=0.01;


#define NEW_FILE_PERM (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH)

enum class SketchType {
    CM = 0,
    COUNT,
    TOWER,
    NITROCM,
    NITROCOUNT
};

class Constraint
{
public:
    double err;
    double prob;
};

class Demand
{
public:
    std::string name;
    SketchType type;
    std::deque<Constraint> constraints;
    double sample_rate;
};

class Config
{
public:
    int type;
    int nrows;
    int ncols;
};

[[maybe_unused]]
static int Open(const char* file, int flag, int perm = NEW_FILE_PERM)
{
    int fd;
    if ((fd=open(file,flag,perm))<0)
    {
        LOG_ERROR("Open file error!\n");
        LOG_ERROR("Can not open file: %s\n", file);
        exit(-1);
    }
    return fd;
}

[[maybe_unused]]
static void Write(int fd, const void* buf, size_t len)
{
    if (write(fd,buf,len)<0)
    {
        LOG_ERROR("Write error!\n");
        exit(-1);
    }
}

inline bool stop(int niter, double prev_avg, double cur_avg)
{
    if (std::fabs(cur_avg-prev_avg) > STOP_THRESHOLD*prev_avg)
        return false;
    else
        return true;
}

inline bool stop(int niter, double prev_avg, double cur_avg, double demand_prob)
{
    if (std::fabs(cur_avg-prev_avg) > STOP_THRESHOLD*prev_avg)
        return false;
    else
        return true;
}

typedef std::chrono::high_resolution_clock::time_point TP;

inline TP now() { return std::chrono::high_resolution_clock::now(); }
