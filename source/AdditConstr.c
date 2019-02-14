//
//  Add_constr.c
//  NumMin
//
//  Created by Alessandro Coretti on 13/02/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "AdditConstr.h"

void ReducedProblem(int N, double x_const, double ** A, double * b, double * x0, int verbose, int werbose, double * xF_AN, double * xF_CG, double * xF_SH, double * xF_BSH) {

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

      b_reduced[i] = 2*x_const*(A_N[i] - A_NN) + b[i] - b_N;
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
