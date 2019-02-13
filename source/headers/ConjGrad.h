//
//  ConjGrad.h
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef ConjGrad_h
#define ConjGrad_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Vector.h"
#include "aux_func.h"
#include "output.h"

double * ConjugateGradient(int N, double ** A, double * b, double * x0, double tol, int maxiter, double * xk);

#endif /* ConjGrad_h */
