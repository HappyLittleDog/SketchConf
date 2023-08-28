#pragma once
#include <deque>
#include "defs.hpp"
#include "streamgen.h"

int poisson_trunc(int l, const Constraint& constraint, StreamGen& stream);