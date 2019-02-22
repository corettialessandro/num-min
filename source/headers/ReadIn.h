//
//  ReadIn.h
//  NumMin
//
//  Created by Alessandro Coretti on 12/13/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef ReadIn_h
#define ReadIn_h

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>

void GetOptions(int argc, char * argv[], int * verbose, int * werbose);
void ReadInput(char * inputfile, int * dimension, int * nblocks, char * Afilename, char * bfilename, char * x0filename, double * tol, int * maxiter, char * mode, double * x_const);

#endif /* ReadIn_h */
