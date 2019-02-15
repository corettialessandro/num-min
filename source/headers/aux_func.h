//
//  aux_func.h
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef aux_func_h
#define aux_func_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Vector.h"

double Norm_inf(int N, double * v);
double Norm_l2(int N, double * v);
double modsq(int N, double * v);
double ScalProd(int N, double * u, double * v);
double * RowbyColProd(int N, double ** M, double * v, double * result);
double * VectorSum(int N, double * u, double * v, double * result);
double * VectorDifference(int N, double * u, double * v, double * result);
double MaxIterError(int N, double * an_sol, double * iter_sol);
double SumComponents(int N, double * x);

#endif /* aux_func_h */
