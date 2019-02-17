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

double ** Reduce_A(int N, double ** A, double ** A_reduced);
double * Reduce_b(int N, double * b, double * A_N, double x_const, double * b_reduced);
double * Reduce_x0(int N, double * x0, double * x0_reduced);
double LastComponent(int N, double * x_others, double x_const);
void CheckAdditionalConstraint(int N, double * x, double x_const, double tol);
void CheckOtherConstraints(int N, double ** A, double * b, double * x, double tol);

#endif /* Reduced_h */
