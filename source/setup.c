//
//  setup.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/17/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "setup.h"

void Initialize(int N, double ** A, double * b, double * x0, char * Afilename, char * bfilename, char * x0filename, int verbose, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH) {

   A = AllocateMatrix(N, N);
   b = AllocateDVector(N);
   x0 = AllocateDVector(N);
   xF_AN = AllocateDVector(N);
   xF_CG = AllocateDVector(N);
   xF_SH = AllocateDVector(N);
   xF_BSH = AllocateDVector(N);

   A = ReadMatrix(N, N, A, Afilename);
   b = ReadDVector(N, b, bfilename);
   x0 = ReadDVector(N, x0, x0filename);

   if (verbose) PrintMatrix(N, N, A, "A");
   if (verbose) PrintVector(N, b, "b");
   if (verbose) PrintVector(N, x0, "x0");

   return;
}

void Finalize(int N, double ** A, double * b, double * x0, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH) {

   FreeMatrix(N, N, A);
   FreeDVector(N, b);
   FreeDVector(N, x0);
   FreeDVector(N, xF_AN);
   FreeDVector(N, xF_CG);
   FreeDVector(N, xF_SH);
   FreeDVector(N, xF_BSH);

   return;
}
