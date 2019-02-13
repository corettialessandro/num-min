//
//  Matrix.h
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef Matrix_h
#define Matrix_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Accelerate/Accelerate.h>

#include "Vector.h"
#include "output.h"

double ** AllocateMatrix(int rows, int columns);
double ** ReadMatrix(int rows, int columns, double ** matrix, char * filepath);
double ** InvertMatrix(int dimension, double ** matrix, double ** inverse);
double * SpectrumMatrix(int dimension, double ** matrix, double ** eigen_matrix, double * eigenvalues);
double ConditionNumber(int dimension, double * eigenvalues);
int CheckDiagonallyDominance(int dimension, double ** matrix, double * diagonal, double * R);
double rcoefficient(int dimension, double ** matrix);
void FreeMatrix(int rows, int columns, double ** matrix);

#endif /* Matrix_h */
