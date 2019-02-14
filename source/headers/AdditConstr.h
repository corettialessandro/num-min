//
//  Add_constr.h
//  NumMin
//
//  Created by Alessandro Coretti on 13/02/19.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef AdditConstr_h
#define AdditConstr_h

#include "Matrix.h"
#include "Vector.h"
#include "ConjGrad.h"
#include "ML_SHAKE.h"
#include "ML_BSHAKE.h"

void ReducedProblem(int N, double x_const, double ** A, double * b, double * x0, int verbose, int werbose, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH);

double ** Reduce_A(int N, double ** A, double ** A_reduced);
double * Reduce_b(int N, double * b, double * A_N, double x_const, double * b_reduced);
double * Reduce_x0(int N, double * x0, double * x0_reduced);

#endif /* AdditConstr_h */
