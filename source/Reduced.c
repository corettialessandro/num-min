//
//  Reduced.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/13/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "Reduced.h"

double ** Reduce_A(int N, double ** A, double ** A_reduced){

   int i, j;
   double A_NN = A[N-1][N-1];
   double * A_N;

   A_N = AllocateDVector(N-1);

   for (i = 0; i < N-1; i++) {

      A_N [i] = A[N-1][i];
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
      printf("[DEBUG] sum = %.12e, x_const = %.12e\n", sum_x, x_const);
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
   }

   discr = Norm_inf(N, sigma);

   if (discr > 1.E1*tol) {
      printf("\nAdditConstr.c -> CheckOtherConstraints() Error: Other constraints not satisfied!\n");
      for (i = 0; i < N; i++) {
         printf("[DEBUG] sigma[%d] = %.12e\n", i, sigma[i]);
      }
      exit(EXIT_FAILURE);
   }

   FreeDVector(N, sigma);

   return;
}
