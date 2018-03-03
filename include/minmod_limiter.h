#ifndef MINMOD_LIMITER_H
#define MINMOD_LIMITER_H

#include <algorithm>
#include <iostream>
#include "para.h"

class minmod_limiter {
  public:
    static void run(SOL& sol, const double h);
    static void run(VEC<EVEC>& sol, const double h);
};

#endif //MINMOD_LIMITER_H

