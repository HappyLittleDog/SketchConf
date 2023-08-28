#pragma once
#include "defs.hpp"
#include "streamgen.h"
using namespace std;

bool satisfy(int ncols, int nrows, const Demand& demand, StreamGen& stream, bool trivial, bool trunc, bool reuse);
pair<int, int> search(StreamGen& stream, const Demand& demand, bool trivial, bool trunc, bool reuse);
pair<int, int> StrawmanSearch(StreamGen& stream, const Demand& demand);
