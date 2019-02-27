//
//  Vector.h
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef Vector_h
#define Vector_h

#include <stdio.h>
#include <stdlib.h>

double * AllocateDVector(int dim);
int * AllocateIVector(int dim);
double * CAllocateDVector(int dim);
double * ReadDVector(int dim, double * vector, char * filepath);
void FreeDVector(int dim, double * vector);
void FreeIVector(int dim, int * vector);

#endif /* Vector_h */
