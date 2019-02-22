//
//  Minimization.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/17/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "Minimization.h"

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

   xF_SH = MasslessShake(N, N, A, b, x0, tol, maxiter, -1, xF_SH);
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

void Constrained(int N, double ** A, double * b, double * x0, int verbose, int werbose, double tol, int maxiter, double x_const, double * xF_SH) {

   xF_SH = MasslessShake(N+1, N, A, b, x0, tol, maxiter, x_const, xF_SH);

   if (verbose) PrintVector(N, xF_SH, "xF_SH");

   return;
}

void Reduced(int N, double x_const, double ** A, double * b, double * x0, int verbose, int werbose, double tol, int maxiter, int nblocks, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH) {

   double red_maxerr_CG, red_maxerr_SH;
   double * red_b, * red_x0;
   double ** red_A, ** red_Ainv;
   clock_t red_an_tstart, red_an_tfinish;
   double red_an_time;

   red_A = AllocateMatrix(N-1, N-1);
   red_b = AllocateDVector(N-1);
   red_x0 = AllocateDVector(N-1);
   red_Ainv = AllocateMatrix(N-1, N-1);

   red_A = Reduce_A(N, A, red_A);
   red_b = Reduce_b(N, b, A[N-1], x_const, red_b);
   red_x0 = Reduce_x0(N, x0, red_x0);

   if (verbose) PrintMatrix(N-1, N-1, red_A, "A Reduced");
   if (verbose) PrintVector(N-1, red_b, "b Reduced");
   if (verbose) PrintVector(N-1, x0, "x0 Reduced");

   if (werbose) WriteMatrix("red_A.out", "output/", N-1, N-1, red_A, "A Reduced");
   if (werbose) WriteVector("red_b.out", "output/", N-1, red_b, "b Reduced");
   if (werbose) WriteVector("red_x0.out", "output/", N-1, red_x0, "x0 Reduced");

   red_an_tstart = clock();
   red_Ainv = InvertMatrix(N-1, red_A, red_Ainv);
   red_an_tfinish = clock();
   red_an_time = (double) (red_an_tfinish - red_an_tstart)/CLOCKS_PER_SEC;
   PrintStats('A', 0, 0, red_an_time, 0, 0);

   xF_AN = RowbyColProd(N-1, red_Ainv, red_b, xF_AN);
   xF_CG = ConjugateGradient(N-1, red_A, red_b, red_x0, tol, maxiter, xF_CG);
   red_maxerr_CG = MaxIterError(N-1, xF_AN, xF_CG);
   xF_SH = MasslessShake(N-1, N-1, red_A, red_b, red_x0, tol, maxiter, -1, xF_SH);
   red_maxerr_SH = MaxIterError(N-1, xF_AN, xF_SH);
   // xF_BSH = MasslessBlockShake(N-1, nblocks, red_A, red_b, red_x0, tol, maxiter, xF_BSH);

   CheckOtherConstraints(N-1, red_A, red_b, xF_AN, tol);
   CheckOtherConstraints(N-1, red_A, red_b, xF_CG, tol);
   CheckOtherConstraints(N-1, red_A, red_b, xF_SH, tol);

   if (verbose) PrintVector(N-1, xF_AN, "xF_AN Reduced");
   if (verbose) PrintVector(N-1, xF_CG, "xF_CG Reduced");
   if (verbose) PrintVector(N-1, xF_SH, "xF_SH Reduced");
   if (werbose) WriteVector("red_xF_AN.out", "output/", N-1, xF_AN, "xF_AN Reduced");
   if (werbose) WriteVector("red_xF_CG.out", "output/", N-1, xF_CG, "xF_CG Reduced");
   if (werbose) WriteVector("red_xF_SH.out", "output/", N-1, xF_SH, "xF_SH Reduced");

   xF_AN[N-1] = LastComponent(N, xF_AN, x_const);
   xF_CG[N-1] = LastComponent(N, xF_CG, x_const);
   xF_SH[N-1] = LastComponent(N, xF_SH, x_const);

   if (verbose) PrintVector(N, xF_AN, "xF_AN");
   if (verbose) PrintVector(N, xF_CG, "xF_CG");
   if (verbose) PrintVector(N, xF_SH, "xF_SH");
   if (werbose) WriteVector("xF_AN.out", "output/", N, xF_AN, "xF_AN");
   if (werbose) WriteVector("xF_CG.out", "output/", N, xF_CG, "xF_CG");
   if (werbose) WriteVector("xF_SH.out", "output/", N, xF_SH, "xF_SH");

   CheckAdditionalConstraint(N, xF_CG, x_const, tol);
   CheckAdditionalConstraint(N, xF_SH, x_const, tol);

   CheckOtherConstraints(N, A, b, xF_CG, tol);
   CheckOtherConstraints(N, A, b, xF_SH, tol);

   FreeMatrix(N-1, N-1, red_A);
   FreeDVector(N-1, red_b);
   FreeDVector(N-1, red_x0);
   FreeMatrix(N-1, N-1, red_Ainv);

   return;
}
