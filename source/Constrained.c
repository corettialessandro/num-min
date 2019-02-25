//
//  Constrained.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/22/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "Constrained.h"

void Constrained(int N, double ** A, double * b, double * x0, int verbose, int werbose, double tol, int maxiter, double x_const, double * xF_AN, double * xF_SH) {

   double maxerr_SH;
   double * constr_b;
   double ** Ainv;
   clock_t constr_an_tstart, constr_an_tfinish;
   double constr_an_time;

   Ainv = AllocateMatrix(N, N);
   constr_b = AllocateDVector(N);

   constr_an_tstart = clock();
   Ainv = InvertMatrix(N, A, Ainv);
   constr_b = Constrain_b(N, Ainv, b, x_const, constr_b);
   xF_AN = RowbyColProd(N, Ainv, constr_b, xF_AN);
   constr_an_tfinish = clock();
   constr_an_time = (double) (constr_an_tfinish - constr_an_tstart)/CLOCKS_PER_SEC;
   PrintStats('A', 0, 0, constr_an_time, 0, 0);
   WriteStats("logfile.out", "output/", 'A', 0, 0, constr_an_time, 0, 0);
   if (verbose) PrintMatrix(N, N, Ainv, "A^{-1}");
   if (verbose) PrintVector(N, constr_b, "b Constrained");
   if (verbose) PrintVector(N, xF_AN, "xF_AN");
   if (werbose) WriteMatrix("Ainv.out", "output/", N, N, Ainv, "A^{-1}");
   if (werbose) WriteVector("constr_b.out", "output/", N, constr_b, "b Constrained");
   if (werbose) WriteVector("xF_AN.out", "output/", N, xF_AN, "xF_AN");

   xF_SH = MasslessShake(N+1, N, A, b, x0, tol, maxiter, x_const, 0, 0, xF_SH);
   maxerr_SH = MaxIterError(N, xF_AN, xF_SH);
   WriteError("logfile.out", "output/", 'S', maxerr_SH);
   if (verbose) PrintVector(N, xF_SH, "xF_SH");
   if (werbose) WriteVector("xF_SH.out", "output/", N, xF_SH, "xF_SH");
   printf("  Maximum error component on iterative solution (SH):\n ");
   printf("  MAX|xf_AN - xf_SH| = %.4e\n\n", maxerr_SH);

   FreeMatrix(N, N, Ainv);
   FreeDVector(N, constr_b);

   return;
}

double * Constrain_b(int N, double ** Ainv, double * b, double x_const, double * constr_b) {

   int i, j;
   double IdAinvb = 0., IdAinvId = 0.;

   for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {

         IdAinvb += Ainv[i][j]*b[j];
         IdAinvId += Ainv[i][j];
      }
   }

   for (i = 0; i < N; i++) {

      constr_b[i] = b[i] + (x_const-IdAinvb)/IdAinvId;
   }

   printf("[DEBUG] gamma_N = %lf\n", (x_const-IdAinvb)/IdAinvId);

   return constr_b;
}
