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

double * Constrain_b(int N, double ** Ainv, double * b, double x_const, double * constr_b);

#endif /* Constrained_h */
