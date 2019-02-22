//
//  Minimization.h
//  NumMin
//
//  Created by Alessandro Coretti on 02/17/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#ifndef Minimization_h
#define Minimization_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "Vector.h"
#include "Matrix.h"
#include "aux_func.h"
#include "output.h"
#include "ConjGrad.h"
#include "ML_SHAKE.h"
#include "ML_BSHAKE.h"
#include "Reduced.h"

void Unconstrained(int N, double ** A, double * b, double * x0, int verbose, int werbose, double tol, int maxiter, int nblocks, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH);
void Constrained();
void Reduced(int N, double x_const, double ** A, double * b, double * x0, int verbose, int werbose, double tol, int maxiter, int nblocks, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH);

#endif /* Minimization_h */
