#ifndef SCALING_LIMITER_H
#define SCALING_LIMITER_H

#include <algorithm>
#include <iostream>
#include "para.h"

#define VECD VEC<double>

class scaling_limiter {
  public:
    //scaling_limiter();
    //scaling_limiter(const VEC<VECD>& av, VEC<VEC<VECD>>& val);
    static void run(const VEC<double>& val, EVEC& sol);
    static void run(const VEC<VEC<VEC<double>>>& val, SOL& sol);
};

#endif //SCALING_LIMITER_H

