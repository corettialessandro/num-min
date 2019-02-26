//
//  ML_SHAKE.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "ML_SHAKE.h"

double * MasslessShake(int N_constr, int N_var, double ** A, double * b, double * x0, double * constr, double tol, int maxiter, double x_const, double * xold) {

    clock_t tstart = clock(), tfinish;
    double time;

    int iter, k, i;
    double discr = 0, gamk;
    double * denom, * sigold;

    denom = AllocateDVector(N_constr);
    sigold = AllocateDVector(N_constr);

    for (i=0; i<N_var; i++) xold[i] = x0[i];

    for (k=0; k<N_constr; k++) {

        denom[k] = Denom(k, N_var, A);
        sigold[k] = Sigma(k, N_var, A, b, xold, constr, x_const);
    }

    discr = Norm_inf(N_constr, sigold);

    for (iter=0; iter<maxiter; iter++) {

        for (k=0; k<N_constr; k++) {

            sigold[k] = Sigma(k, N_var, A, b, xold, constr, x_const);

            if (fabs(sigold[k]) > tol) {

                if (discr < fabs(sigold[k])) discr = fabs(sigold[k]);

                gamk = sigold[k]/denom[k];

                if (k >= N_var) {

                   constr[k - N_var] -= gamk;
                }

                for (i=0; i<N_var; i++) {

                    if (k < N_var) xold[i] -= gamk*A[k][i];
                }
            }
        }

        discr = Norm_inf(N_constr, sigold);

        WriteConvergence('S', "SH_convergence.out", "output/", iter, discr, maxiter, tol);

        if (discr < tol) break;
    }

    // for (k = 0; k < N_constr; k++) printf("sigma[%d] = %.12e\n", k, sigold[k] + constr[k]);

    if (iter == maxiter){

        printf("\nML_SHAKE.c -> MasslessShake() Error: Maximum number of iterations reached!\n");
//        exit(EXIT_FAILURE);
    }

    FreeDVector(N_constr, denom);
    FreeDVector(N_constr, sigold);

    tfinish = clock();

    time = (double) (tfinish - tstart)/CLOCKS_PER_SEC;

    PrintStats('S', 0, iter, time, discr, tol);
    WriteStats("logfile.out", "output/", 'S', 0, iter, time, discr, tol);

    return xold;
}

double Sigma(int k, int N_var, double ** A, double * b, double * x, double * constr, double x_const) {

   double sigma_k;

   int i;

   if (k < N_var) {

      sigma_k = -b[k] - constr[k];
      for (i = 0; i < N_var; i++) sigma_k += A[k][i] * x[i];

   } else {

      sigma_k = -x_const;
      for (i = 0; i < N_var; i++) sigma_k += x[i];
   }

   return sigma_k;
}

double Denom(int k, int N_var, double ** A) {

   double denom_k;

   int i;

   if (k < N_var) {

      denom_k = 0.;
      for (i = 0; i < N_var; i++) denom_k += A[k][i]*A[k][i];

   } else {

      denom_k = N_var;
   }

   return denom_k;
}
