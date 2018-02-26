/**
 * @file BasFun.h
 * @brief Basis function on [-1,1], i.e. Legendre polynomial not normalized.
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-04
 */

#ifndef BASFUN_H
#define BASFUN_H

#include <vector>
#include <iostream>
#include "para.h"

VEC<double> Poly(double x);
VEC<double> PolyG(double x);

#endif //BASFUN_H

