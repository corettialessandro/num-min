//
//  Tensor.h
//  NumMin
//
//  Created by Alessandro Coretti on 12/20/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef Tensor_h
#define Tensor_h

#include <stdio.h>
#include <stdlib.h>

double *** AllocateTensor(int rows, int columns, int sections);
void FreeTensor(int rows, int columns, int sections, double *** tensor);

#endif /* Tensor_h */
