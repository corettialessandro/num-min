//
//  Add_constr.c
//  NumMin
//
//  Created by Alessandro Coretti on 13/02/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "AdditConstr.h"

void ReducedProblem(int N, double x_const, double ** A, double * b, double * x0, int verbose, int werbose, double tol, int maxiter, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH) {

   double * red_b, * red_x0, * red_xF_AN, * red_xF_CG, * red_xF_SH, * red_xF_BSH;
   double ** red_A, ** red_Ainv;

   red_A = AllocateMatrix(N-1, N-1);
   red_b = AllocateDVector(N-1);
   red_x0 = AllocateDVector(N-1);
   red_Ainv = AllocateMatrix(N-1, N-1);
   red_xF_AN = AllocateDVector(N-1);
   red_xF_CG = AllocateDVector(N-1);
   red_xF_SH = AllocateDVector(N-1);
   red_xF_BSH = AllocateDVector(N-1);

   red_A = Reduce_A(N, A, red_A);
   red_b = Reduce_b(N, b, A[N-1], x_const, red_b);
   red_x0 = Reduce_x0(N, x0, red_x0);

   if (verbose) PrintMatrix(N-1, N-1, red_A, "A Reduced");
   if (verbose) PrintVector(N-1, red_b, "b Reduced");
   if (verbose) PrintVector(N-1, x0, "x0 Reduced");

   if (werbose) WriteMatrix("red_A.out", "output/", N-1, N-1, red_A, "A Reduced");
   if (werbose) WriteVector("red_b.out", "output/", N-1, red_b, "b Reduced");
   if (werbose) WriteVector("red_x0.out", "output/", N-1, x0, "x0 Reduced");

   xF_CG = ConjugateGradient(N-1, red_A, red_b, red_x0, tol, maxiter, xF_CG);
   xF_SH = MasslessShake(N-1, red_A, red_b, red_x0, tol, maxiter, xF_SH);
   // xF_BSH = MasslessBlockShake(N-1, nblocks, red_A, red_b, red_x0, tol, maxiter, xF_BSH);

   CheckOtherConstraints(N-1, red_A, red_b, xF_CG, tol);
   CheckOtherConstraints(N-1, red_A, red_b, xF_SH, tol);

   xF_CG[N-1] = LastComponent(N, xF_CG, x_const);
   xF_SH[N-1] = LastComponent(N, xF_SH, x_const);

   if (verbose) PrintVector(N, xF_CG, "xF_CG");
   if (verbose) PrintVector(N, xF_SH, "xF_SH");
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
   FreeDVector(N-1, red_xF_AN);
   FreeDVector(N-1, red_xF_CG);
   FreeDVector(N-1, red_xF_SH);
   FreeDVector(N-1, red_xF_BSH);

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

         A_reduced[i][j] = A_NN - 2*A_N[j] + A[i][j];
      }
   }

   FreeDVector(N-1, A_N);

   return A_reduced;
}

double * Reduce_b(int N, double * b, double * A_N, double x_const, double * b_reduced){

   int i;
   double b_N = b[N-1], A_NN = A_N[N-1];

   for (i = 0; i < N-1; i++) {

      b_reduced[i] = -(2*x_const*(A_N[i] - A_NN) - b[i] + b_N);
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

double LastComponent(int N, double * x_others, double x_const){

   int i;
   double xlast = 0;

   for (i = 0; i < N-1; i++) {

      xlast -= x_others[i];
   }

   xlast += x_const;

   return xlast;
}

void CheckAdditionalConstraint(int N, double * x, double x_const, double tol) {

   int i;
   double sum_x = 0;

   for (i = 0; i < N; i++) {

      sum_x += x[i];
   }

   if (fabs(sum_x - x_const) > tol) {

      printf("\nAdditConstr.c -> CheckAdditionalConstraint() Error: Additional constraint not satisfied!\n");
      exit(EXIT_FAILURE);
   }

   return;
}

void CheckOtherConstraints(int N, double ** A, double * b, double * x, double tol) {

   int k, i;
   double discr = 0.;
   double * sigma;

   sigma = AllocateDVector(N);

   for (k = 0; k < N; k++) {

      sigma[k] = -b[k];

      for (i = 0; i < N; i++) {

         sigma[k] += A[k][i]*x[i];
      }

      printf("%lf\n", sigma[k]);
   }

   discr = Norm_inf(N, sigma);

   if (discr > tol) {

      printf("\nAdditConstr.c -> CheckOtherConstraints() Error: Other constraints not satisfied!\n");
      exit(EXIT_FAILURE);
   }

   FreeDVector(N, sigma);

   return;
}
