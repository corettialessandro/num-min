//
//  Unconstrained.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/17/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "Unconstrained.h"

void Unconstrained(int N_var, int N_constr, double ** A, double * b, double * x0, double * constr, int verbose, int werbose, double tol, int maxiter, int nblocks, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH) {

   double maxerr_CG, maxerr_SH, maxerr_BSH;
   double lambda = 0;
   double * bconstr;
   double ** Ainv;
   clock_t an_tstart, an_tfinish;
   double an_time;

   Ainv = AllocateMatrix(N_var, N_var);
   bconstr = AllocateDVector(N_var);

   an_tstart = clock();
   Ainv = InvertMatrix(N_var, A, Ainv);
   if (N_constr - N_var == 1) lambda = Lambda(N_var, Ainv, b, constr[N_constr-1]);
   bconstr = VectorScalarSum(N_var, b, lambda, bconstr);
   xF_AN = RowbyColProd(N_var, Ainv, bconstr, xF_AN);
   an_tfinish = clock();
   an_time = (double) (an_tfinish - an_tstart)/CLOCKS_PER_SEC;
   PrintStats('A', 0, 0, an_time, 0, 0);
   WriteStats("logfile.out", "output/", 'A', 0, 0, an_time, 0, 0);
   if (verbose) PrintMatrix(N_var, N_var, Ainv, "A^{-1}");
   if (verbose) PrintVector(N_var, xF_AN, "xF_AN");
   if (werbose) WriteMatrix("Ainv.out", "output/", N_var, N_var, Ainv, "A^{-1}");
   if (werbose) WriteVector("xF_AN.out", "output/", N_var, xF_AN, "xF_AN");

   CheckConstraints(N_var, N_constr, A, b, xF_AN, constr);

   // xF_CG = ConjugateGradient(N_var, A, b, x0, tol, maxiter, xF_CG);
   // maxerr_CG = MaxIterError(N_var, xF_AN, xF_CG);
   // WriteError("logfile.out", "output/", 'C', maxerr_CG);
   // if (verbose) PrintVector(N_var, xF_CG, "xF_CG");
   // if (werbose) WriteVector("xF_CG.out", "output/", N_var, xF_CG, "xF_CG");
   // printf("  Maximum error component on iterative solution (CG):\n");
   // printf("  MAX|xf_AN - xf_CG| = %.4e\n\n", maxerr_CG);

   xF_SH = MasslessShake(N_var, N_constr, A, b, x0, constr, tol, maxiter, xF_SH);
   maxerr_SH = MaxIterError(N_var, xF_AN, xF_SH);
   WriteError("logfile.out", "output/", 'S', maxerr_SH);
   if (verbose) PrintVector(N_var, xF_SH, "xF_SH");
   if (werbose) WriteVector("xF_SH.out", "output/", N_var, xF_SH, "xF_SH");
   printf("  Maximum error component on iterative solution (SH):\n ");
   printf("  MAX|xf_AN - xf_SH| = %.4e\n\n", maxerr_SH);

   if (verbose) CheckConstraints(N_var, N_constr, A, b, xF_SH, constr);

   // xF_BSH = MasslessBlockShake(N_var, nblocks, A, b, x0, tol, maxiter, xF_BSH);
   // maxerr_BSH = MaxIterError(N_var, xF_AN, xF_BSH);
   // WriteError("logfile.out", "output/", 'B', maxerr_BSH);
   // if (verbose) PrintVector(N_var, xF_BSH, "xF_BSH");
   // if (werbose) WriteVector("xF_BSH.out", "output/", N_var, xF_BSH, "xF_BSH");
   // printf("  Maximum error component on iterative solution (BSH):\n");
   // printf("  MAX|xf_AN - xf_BSH| = %.4e\n\n", maxerr_BSH);

   FreeMatrix(N_var, N_var, Ainv);
   FreeDVector(N_var, bconstr);

   return;
}

double Lambda(int N, double ** Ainv, double * b, double x_const) {

   int i, j;
   double IdAinvb = 0., IdAinvId = 0.;

   for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {

         IdAinvb += Ainv[i][j]*b[j];
         IdAinvId += Ainv[i][j];
      }
   }

   return (x_const-IdAinvb)/IdAinvId;
}
