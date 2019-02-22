//
//  Constrained.h
//  NumMin
//
//  Created by Alessandro Coretti on 02/22/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#ifndef Constrained_h
#define Constrained_h

#include <stdio.h>

#include "Matrix.h"
#include "Vector.h"
#include "aux_func.h"
#include "output.h"
#include "ConjGrad.h"
#include "ML_SHAKE.h"
#include "ML_BSHAKE.h"

void Constrained(int N, double ** A, double * b, double * x0, int verbose, int werbose, double tol, int maxiter, double x_const, double * xF_AN, double * xF_SH);
double * Constrain_b(int N, double ** Ainv, double * b, double x_const, double * constr_b);

#endif /* Constrained_h */
