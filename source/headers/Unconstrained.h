//
//  Unconstrained.h
//  NumMin
//
//  Created by Alessandro Coretti on 02/17/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#ifndef Unconstrained_h
#define Unconstrained_h

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

void Unconstrained(int N_var, int N_constr, double ** A, double * b, double * x0, double * constr, int verbose, int werbose, double tol, int maxiter, int nblocks, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH);

#endif /* Unconstrained_h */
