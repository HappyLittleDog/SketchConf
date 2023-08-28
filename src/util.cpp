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

int poisson_trunc(int l, const Constraint& constraint, StreamGen& stream)
{
    double nflows=stream.getTotalFlows();
    double lambda=nflows/l;
    double rst=5*sqrt(lambda);

    double xi = sqrt(lambda);
    xi = 5.0 / (5.0 + xi);

    double K = -log(1-xi)/xi;
    double delta = 1.0 / (constraint.prob * RELATIVE_ERROR);
    double rst2 = log(delta) / (K-1);

    return lambda + max(rst, rst2); 
}