//
//  output.h
//  NumMin
//
//  Created by Alessandro Coretti on 12/13/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef output_h
#define output_h

#include <stdio.h>
#include <stdlib.h>

#include "Vector.h"
#include "aux_func.h"

void PrintSetup(char * inputfile, int dimension, int nblocks, char * Afilename, char * bfilename, char * x0filename, double tol, int maxiter, int additional_constraint, double x_const);
void PrintMatrix(int rows, int columns, double ** matrix, char * label);
void PrintVector(int dim, double * vector, char * label);
void PrintStats(char method, int nblocks, int iter, double exectime, double discr, double tol);

void WriteSetup(char * outputfile, char * outputpath, int dimension, int nblocks, char * Afilename, char * bfilename, char * x0filename, double tol, int maxiter, int additional_constraint, double x_const);
void WriteMatrix(char * outputfile, char * outputpath, int rows, int columns, double ** matrix, char * label);
void WriteVector(char * outputfile, char * outputpath, int dim, double * vector, char * label);
void WriteQFProp(char * outputfile, char * outputpath, int dd, double r, double scn);
void WriteStats(char * outputfile, char * outputpath, char method, int nblocks, int iter, double exectime, double discr, double tol);
void WriteError(char * outputfile, char * outputpath, char method, double maxerr);
void WriteConvergence(char method, char * outputfile, char * outputpath, int iter, double discr, int maxiter, double tol);

#endif /* output_h */
