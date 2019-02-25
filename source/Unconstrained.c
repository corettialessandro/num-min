//
//  Unconstrained.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/17/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "Unconstrained.h"

void Unconstrained(int N, double ** A, double * b, double * x0, int verbose, int werbose, double tol, int maxiter, int nblocks, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH) {

   double maxerr_CG, maxerr_SH, maxerr_BSH;
   double ** Ainv;
   clock_t an_tstart, an_tfinish;
   double an_time;

   Ainv = AllocateMatrix(N, N);

   an_tstart = clock();
   Ainv = InvertMatrix(N, A, Ainv);
   xF_AN = RowbyColProd(N, Ainv, b, xF_AN);
   an_tfinish = clock();
   an_time = (double) (an_tfinish - an_tstart)/CLOCKS_PER_SEC;
   PrintStats('A', 0, 0, an_time, 0, 0);
   WriteStats("logfile.out", "output/", 'A', 0, 0, an_time, 0, 0);
   if (verbose) PrintMatrix(N, N, Ainv, "A^{-1}");
   if (verbose) PrintVector(N, xF_AN, "xF_AN");
   if (werbose) WriteMatrix("Ainv.out", "output/", N, N, Ainv, "A^{-1}");
   if (werbose) WriteVector("xF_AN.out", "output/", N, xF_AN, "xF_AN");

   xF_CG = ConjugateGradient(N, A, b, x0, tol, maxiter, xF_CG);
   maxerr_CG = MaxIterError(N, xF_AN, xF_CG);
   WriteError("logfile.out", "output/", 'C', maxerr_CG);
   if (verbose) PrintVector(N, xF_CG, "xF_CG");
   if (werbose) WriteVector("xF_CG.out", "output/", N, xF_CG, "xF_CG");
   printf("  Maximum error component on iterative solution (CG):\n");
   printf("  MAX|xf_AN - xf_CG| = %.4e\n\n", maxerr_CG);

   xF_SH = MasslessShake(N, N, A, b, x0, tol, maxiter, 0, +1, -1, xF_SH);
   maxerr_SH = MaxIterError(N, xF_AN, xF_SH);
   WriteError("logfile.out", "output/", 'S', maxerr_SH);
   if (verbose) PrintVector(N, xF_SH, "xF_SH");
   if (werbose) WriteVector("xF_SH.out", "output/", N, xF_SH, "xF_SH");
   printf("  Maximum error component on iterative solution (SH):\n ");
   printf("  MAX|xf_AN - xf_SH| = %.4e\n\n", maxerr_SH);

   xF_BSH = MasslessBlockShake(N, nblocks, A, b, x0, tol, maxiter, xF_BSH);
   maxerr_BSH = MaxIterError(N, xF_AN, xF_BSH);
   WriteError("logfile.out", "output/", 'B', maxerr_BSH);
   if (verbose) PrintVector(N, xF_BSH, "xF_BSH");
   if (werbose) WriteVector("xF_BSH.out", "output/", N, xF_BSH, "xF_BSH");
   printf("  Maximum error component on iterative solution (BSH):\n");
   printf("  MAX|xf_AN - xf_BSH| = %.4e\n\n", maxerr_BSH);

   if (verbose) printf("  The sum of the components of the minimum point is %.12lf\n\n", SumComponents(N, xF_AN));

   FreeMatrix(N, N, Ainv);

   return;
}
