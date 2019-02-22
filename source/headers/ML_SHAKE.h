//
//  ML_SHAKE.h
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef ML_SHAKE_h
#define ML_SHAKE_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "Matrix.h"
#include "Vector.h"
#include "aux_func.h"
#include "output.h"

double * MasslessShake(int N_constr, int N_var, double ** A, double * b, double * x0, double tol, int maxiter, double x_const, double * xold);
double Sigma(int k, int N_var, double ** A, double * b, double * x, double x_const, double gamma_N);
double Denom(int k, int N_var, double ** A);

#endif /* ML_SHAKE_h */
