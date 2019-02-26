//
//  Reduced.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/13/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "Reduced.h"

void Reduced(int N, double x_const, double ** A, double * b, double * x0, double * constr, int verbose, int werbose, double tol, int maxiter, int nblocks, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH) {

   double red_maxerr_CG, red_maxerr_SH;
   double maxerr_CG, maxerr_SH;
   double * red_b, * red_x0, * red_constr;
   double ** red_A, ** red_Ainv;
   clock_t red_an_tstart, red_an_tfinish;
   double red_an_time;

   red_A = AllocateMatrix(N-1, N-1);
   red_b = AllocateDVector(N-1);
   red_x0 = AllocateDVector(N-1);
   red_constr = AllocateDVector(N-1);
   red_Ainv = AllocateMatrix(N-1, N-1);

   red_A = Reduce_A(N, A, red_A);
   red_b = Reduce_b(N, b, A[N-1], x_const, red_b);
   red_x0 = Reduce_x0(N, x0, red_x0);
   red_constr = Reduce_constr(N, constr, red_constr);

   if (verbose) PrintMatrix(N-1, N-1, red_A, "A Reduced");
   if (verbose) PrintVector(N-1, red_b, "b Reduced");
   if (verbose) PrintVector(N-1, red_x0, "x0 Reduced");
   if (verbose) PrintVector(N-1, red_constr, "constr Reduced");

   if (werbose) WriteMatrix("red_A.out", "output/", N-1, N-1, red_A, "A Reduced");
   if (werbose) WriteVector("red_b.out", "output/", N-1, red_b, "b Reduced");
   if (werbose) WriteVector("red_x0.out", "output/", N-1, red_x0, "x0 Reduced");
   if (werbose) WriteVector("red_constr.out", "output/", N-1, red_constr, "constr Reduced");

   red_an_tstart = clock();
   red_Ainv = InvertMatrix(N-1, red_A, red_Ainv);
   red_an_tfinish = clock();
   red_an_time = (double) (red_an_tfinish - red_an_tstart)/CLOCKS_PER_SEC;
   PrintStats('A', 0, 0, red_an_time, 0, 0);

   xF_AN = RowbyColProd(N-1, red_Ainv, red_b, xF_AN);
   xF_CG = ConjugateGradient(N-1, red_A, red_b, red_x0, tol, maxiter, xF_CG);
   red_maxerr_CG = MaxIterError(N-1, xF_AN, xF_CG);
   xF_SH = MasslessShake(N-1, N-1, red_A, red_b, red_x0, red_constr, tol, maxiter, xF_SH);
   red_maxerr_SH = MaxIterError(N-1, xF_AN, xF_SH);

   // if (verbose) CheckOtherConstraints(N-1, red_A, red_b, xF_AN);
   // if (verbose) CheckOtherConstraints(N-1, red_A, red_b, xF_CG);
   // if (verbose) CheckOtherConstraints(N-1, red_A, red_b, xF_SH);

   if (verbose) PrintVector(N-1, xF_AN, "xF_AN Reduced");
   if (verbose) PrintVector(N-1, xF_CG, "xF_CG Reduced");
   if (verbose) PrintVector(N-1, xF_SH, "xF_SH Reduced");
   if (werbose) WriteVector("red_xF_AN.out", "output/", N-1, xF_AN, "xF_AN Reduced");
   if (werbose) WriteVector("red_xF_CG.out", "output/", N-1, xF_CG, "xF_CG Reduced");
   if (werbose) WriteVector("red_xF_SH.out", "output/", N-1, xF_SH, "xF_SH Reduced");

   xF_AN[N-1] = LastComponent(N, xF_AN, x_const);

   xF_CG[N-1] = LastComponent(N, xF_CG, x_const);
   maxerr_CG = MaxIterError(N, xF_AN, xF_CG);
   printf("  Maximum error component on iterative solution (CG):\n");
   printf("  MAX|xf_AN - xf_CG| = %.4e\n\n", maxerr_CG);

   xF_SH[N-1] = LastComponent(N, xF_SH, x_const);
   maxerr_SH = MaxIterError(N, xF_AN, xF_SH);
   printf("  Maximum error component on iterative solution (SH):\n ");
   printf("  MAX|xf_AN - xf_SH| = %.4e\n\n", maxerr_SH);

   if (verbose) PrintVector(N, xF_AN, "xF_AN");
   if (verbose) PrintVector(N, xF_CG, "xF_CG");
   if (verbose) PrintVector(N, xF_SH, "xF_SH");
   if (werbose) WriteVector("xF_AN.out", "output/", N, xF_AN, "xF_AN");
   if (werbose) WriteVector("xF_CG.out", "output/", N, xF_CG, "xF_CG");
   if (werbose) WriteVector("xF_SH.out", "output/", N, xF_SH, "xF_SH");

   // if (verbose) CheckAdditionalConstraint(N, xF_CG, constr[N_constr - 1]);
   // if (verbose) CheckAdditionalConstraint(N, xF_SH, constr[N_constr - 1]);
   //
   // if (verbose) CheckOtherConstraints(N, A, b, xF_CG);
   // if (verbose) CheckOtherConstraints(N, A, b, xF_SH);

   FreeMatrix(N-1, N-1, red_A);
   FreeDVector(N-1, red_b);
   FreeDVector(N-1, red_x0);
   FreeMatrix(N-1, N-1, red_Ainv);

   return;
}

double ** Reduce_A(int N, double ** A, double ** A_reduced){

   int i, j;
   double A_NN = A[N-1][N-1];
   double * A_N;

   A_N = AllocateDVector(N-1);

   for (i = 0; i < N-1; i++) {

      A_N[i] = A[N-1][i];
   }

   for (i = 0; i < N-1; i++) {
      for (j = 0; j < N-1; j++) {

         A_reduced[i][j] = A_NN - A_N[j] - A_N[i] + A[i][j];
      }
   }

   FreeDVector(N-1, A_N);

   return A_reduced;
}

double * Reduce_b(int N, double * b, double * A_N, double x_const, double * b_reduced){

   int i;
   double b_N = b[N-1], A_NN = A_N[N-1];

   for (i = 0; i < N-1; i++) {

      b_reduced[i] = x_const*(A_NN - A_N[i]) - b_N + b[i];
   }

   return b_reduced;
}

double * Reduce_x0(int N, double * x0, double * x0_reduced){

   int i;

   for (i = 0; i < N-1; i++) {

      x0_reduced[i] = x0[i];
   }

   return x0_reduced;
}

double * Reduce_constr(int N, double * constr, double * constr_reduced){

   int i;

   for (i = 0; i < N-1; i++) {

      constr_reduced[i] = 0.;
   }

   return constr_reduced;
}

double LastComponent(int N, double * x_others, double x_const){

   int i;
   double xlast = 0;

   for (i = 0; i < N-1; i++) {

      xlast -= x_others[i];
   }

   xlast += x_const;

   return xlast;
}
