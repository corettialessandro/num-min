//
//  setup.h
//  NumMin
//
//  Created by Alessandro Coretti on 02/17/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#ifndef setup_h
#define setup_h

#include <stdio.h>

#include <Matrix.h>
#include <Vector.h>
#include <output.h>

void Initialize(int N, double ** A, double * b, double * x0, char * Afilename, char * bfilename, char * x0filename, int verbose, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH);
void Finalize(int N, double ** A, double * b, double * x0, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH);

#endif /* setup_h */
