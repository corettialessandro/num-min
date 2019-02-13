//
//  ConjGrad.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "ConjGrad.h"

double * ConjugateGradient(int N, double ** A, double * b, double * x0, double tol, int maxiter, double * xk) {
    
    clock_t tstart = clock(), tfinish;
    double time;
    
    int i, j, iter = 0;
    double alphak, betakp1, discr, rkdotrk, rkp1dotrkp1;
    double * xkp1, * rk, * rkp1, * dk, * dkp1, * Adk;

    xkp1 = AllocateDVector(N);
    rk = AllocateDVector(N);
    rkp1 = AllocateDVector(N);
    dk = AllocateDVector(N);
    dkp1 = AllocateDVector(N);
    Adk = AllocateDVector(N);

    for (i=0; i<N; i++) {
        
        xk[i] = x0[i];
    }
    
    for (i=0; i<N; i++) {

        rk[i] = -b[i];
        
        for (j=0; j<N; j++) {
            
            rk[i] += A[i][j]*xk[j];
        }
        
        dk[i] = -rk[i];
    }
    
    discr = Norm_inf(N, rk);
    rkdotrk = modsq(N, rk);
    
    WriteConvergence('C', "CG_convergence.out", "output/", iter, discr, maxiter, tol);

    while (discr > tol) {
        
        iter ++;
        
        if (iter == maxiter){
            printf("\nConjGrad.c -> ConjugateGradient() Error: Maximum number of iterations reached!\n");
            exit(EXIT_FAILURE);
        }

        for (i=0; i<N; i++) {
            
            Adk[i] = 0.;

            for (j=0; j<N; j++) {
                
                Adk[i] += A[i][j]*dk[j];
            }
        }
        
        alphak = modsq(N, rk)/ScalProd(N, dk, Adk);
        
        for (i=0; i<N; i++) {
            
            xkp1[i] = xk[i] + alphak*dk[i];
            rkp1[i] = rk[i] + alphak*Adk[i];
        }
        
        if (iter%10 == 0) {

            for (i=0; i<N; i++) {

                rkp1[i] = -b[i];

                for (j=0; j<N; j++) {

                    rkp1[i] += A[i][j]*xkp1[j];
                }
            }
        }
        
        rkp1dotrkp1 = modsq(N, rkp1);
        betakp1 = rkp1dotrkp1/rkdotrk;
        
        for (i=0; i<N; i++) {
            
            dkp1[i] = -rkp1[i] + betakp1*dk[i];
        }
        
        rkdotrk = rkp1dotrkp1;
        for (i=0; i<N; i++) {
            
            xk[i] = xkp1[i];
            rk[i] = rkp1[i];
            dk[i] = dkp1[i];
        }
        
        discr = Norm_inf(N, rk);
        
        WriteConvergence('C', "CG_convergence.out", "output/", iter, discr, maxiter, tol);
    }
    
    FreeDVector(N, xkp1);
    FreeDVector(N, rk);
    FreeDVector(N, rkp1);
    FreeDVector(N, dk);
    FreeDVector(N, dkp1);
    FreeDVector(N, Adk);

    tfinish = clock();
    
    time = (double) (tfinish - tstart)/CLOCKS_PER_SEC;
    
    PrintStats('C', 0, iter, time, discr, tol);
    WriteStats("logfile.out", "output/", 'C', 0, iter, time, discr, tol);
    
    return xk;
}
