//
//  ML_BSHAKE.h
//  NumMin
//
//  Created by Alessandro Coretti on 12/21/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef ML_BSHAKE_h
#define ML_BSHAKE_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "Matrix.h"
#include "Vector.h"
#include "Tensor.h"
#include "aux_func.h"
#include "output.h"

double * MasslessBlockShake(int N, int nblocks, double ** A, double * b, double * x0, double tol, int maxiter, double * xold);

#endif /* ML_BSHAKE_h */
