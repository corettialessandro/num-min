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
   double ** Ainv;
   clock_t an_tstart, an_tfinish;
   double an_time;

   Ainv = AllocateMatrix(N_var, N_var);

   an_tstart = clock();
   Ainv = InvertMatrix(N_var, A, Ainv);
   xF_AN = RowbyColProd(N_var, Ainv, b, xF_AN);
   an_tfinish = clock();
   an_time = (double) (an_tfinish - an_tstart)/CLOCKS_PER_SEC;
   PrintStats('A', 0, 0, an_time, 0, 0);
   WriteStats("logfile.out", "output/", 'A', 0, 0, an_time, 0, 0);
   if (verbose) PrintMatrix(N_var, N_var, Ainv, "A^{-1}");
   if (verbose) PrintVector(N_var, xF_AN, "xF_AN");
   if (werbose) WriteMatrix("Ainv.out", "output/", N_var, N_var, Ainv, "A^{-1}");
   if (werbose) WriteVector("xF_AN.out", "output/", N_var, xF_AN, "xF_AN");

   xF_CG = ConjugateGradient(N_var, A, b, x0, tol, maxiter, xF_CG);
   maxerr_CG = MaxIterError(N_var, xF_AN, xF_CG);
   WriteError("logfile.out", "output/", 'C', maxerr_CG);
   if (verbose) PrintVector(N_var, xF_CG, "xF_CG");
   if (werbose) WriteVector("xF_CG.out", "output/", N_var, xF_CG, "xF_CG");
   printf("  Maximum error component on iterative solution (CG):\n");
   printf("  MAX|xf_AN - xf_CG| = %.4e\n\n", maxerr_CG);

   xF_SH = MasslessShake(N_constr, N_var, A, b, x0, constr, tol, maxiter, 0, xF_SH);
   maxerr_SH = MaxIterError(N_var, xF_AN, xF_SH);
   WriteError("logfile.out", "output/", 'S', maxerr_SH);
   if (verbose) PrintVector(N_var, xF_SH, "xF_SH");
   if (werbose) WriteVector("xF_SH.out", "output/", N_var, xF_SH, "xF_SH");
   printf("  Maximum error component on iterative solution (SH):\n ");
   printf("  MAX|xf_AN - xf_SH| = %.4e\n\n", maxerr_SH);

   xF_BSH = MasslessBlockShake(N_var, nblocks, A, b, x0, tol, maxiter, xF_BSH);
   maxerr_BSH = MaxIterError(N_var, xF_AN, xF_BSH);
   WriteError("logfile.out", "output/", 'B', maxerr_BSH);
   if (verbose) PrintVector(N_var, xF_BSH, "xF_BSH");
   if (werbose) WriteVector("xF_BSH.out", "output/", N_var, xF_BSH, "xF_BSH");
   printf("  Maximum error component on iterative solution (BSH):\n");
   printf("  MAX|xf_AN - xf_BSH| = %.4e\n\n", maxerr_BSH);

   if (verbose) printf("  The sum of the components of the minimum point is %.12lf\n\n", SumComponents(N_var, xF_AN));

   FreeMatrix(N_var, N_var, Ainv);

   return;
}
