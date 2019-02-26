//
//  Reduced.h
//  NumMin
//
//  Created by Alessandro Coretti on 02/13/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#ifndef Reduced_h
#define Reduced_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "aux_func.h"
#include "Matrix.h"
#include "Vector.h"
#include "ConjGrad.h"
#include "ML_SHAKE.h"
#include "ML_BSHAKE.h"

void Reduced(int N, double x_const, double ** A, double * b, double * x0, double * constr, int verbose, int werbose, double tol, int maxiter, int nblocks, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH);
double ** Reduce_A(int N, double ** A, double ** A_reduced);
double * Reduce_b(int N, double * b, double * A_N, double x_const, double * b_reduced);
double * Reduce_x0(int N, double * x0, double * x0_reduced);
double * Reduce_constr(int N, double * constr, double * constr_reduced);
double LastComponent(int N, double * x_others, double x_const);

#endif /* Reduced_h */
