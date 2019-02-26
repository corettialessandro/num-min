//
//  setup.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/17/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "setup.h"

void Initialize(int N_var, int N_constr, double *** A, double ** b, double ** x0, double ** constr, char * Afilename, char * bfilename, char * x0filename, char * constrfilename, int verbose, double ** xF_AN, double ** xF_CG, double ** xF_SH, double ** xF_BSH) {

   *A = AllocateMatrix(N_var, N_var);
   *b = AllocateDVector(N_var);
   *x0 = AllocateDVector(N_var);
   *constr = AllocateDVector(N_constr);
   *xF_AN = AllocateDVector(N_var);
   *xF_CG = AllocateDVector(N_var);
   *xF_SH = AllocateDVector(N_var);
   *xF_BSH = AllocateDVector(N_var);

   *A = ReadMatrix(N_var, N_var, *A, Afilename);
   *b = ReadDVector(N_var, *b, bfilename);
   *x0 = ReadDVector(N_var, *x0, x0filename);
   *constr = ReadDVector(N_constr, *constr, constrfilename);

   if (verbose) PrintMatrix(N_var, N_var, *A, "A");
   if (verbose) PrintVector(N_var, *b, "b");
   if (verbose) PrintVector(N_var, *x0, "x0");
   if (verbose) PrintVector(N_constr, *constr, "constr");

   return;
}

void Finalize(int N_var, int N_constr, double *** A, double ** b, double ** x0, double ** constr, double ** xF_AN, double ** xF_CG, double ** xF_SH, double ** xF_BSH) {

   FreeMatrix(N_var, N_var, *A);
   FreeDVector(N_var, *b);
   FreeDVector(N_var, *x0);
   FreeDVector(N_constr, *constr);
   FreeDVector(N_var, *xF_AN);
   FreeDVector(N_var, *xF_CG);
   FreeDVector(N_var, *xF_SH);
   FreeDVector(N_var, *xF_BSH);

   return;
}
