//
//  ML_SHAKE.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "ML_SHAKE.h"

double * MasslessShake(int N, double ** A, double * b, double * x0, double tol, int maxiter, double * xold) {

    clock_t tstart = clock(), tfinish;
    double time;
    
    int iter, k, i;
    double discr = 0, gamk;
    double * denom, * sigold;
    
    denom = AllocateDVector(N);
    sigold = AllocateDVector(N);
    
    for (i=0; i<N; i++) xold[i] = x0[i];
    
    for (k=0; k<N; k++) {
        
        denom[k] = 0.;
        sigold[k] = -b[k];
        
        for (i=0; i<N; i++) {
            
            denom[k] += A[k][i]*A[k][i];
            sigold[k] += A[k][i]*xold[i];
        }
    }

    discr = Norm_inf(N, sigold);

    for (iter=0; iter<maxiter; iter++) {
        
        for (k=0; k<N; k++) {
            
            sigold[k] = -b[k];
            
            for (i=0; i<N; i++) sigold[k] += A[k][i]*xold[i];
            
            if (fabs(sigold[k]) > tol) {
                
                if (discr < fabs(sigold[k])) discr = fabs(sigold[k]);
                
                gamk = sigold[k]/denom[k];
                
                for (i=0; i<N; i++) {
                    
                    xold[i] -= gamk*A[k][i];
                }
            }
        }
        
        discr = Norm_inf(N, sigold);

        WriteConvergence('S', "SH_convergence.out", "output/", iter, discr, maxiter, tol);
        
        if (discr < tol) break;
    }
    
    if (iter == maxiter){
        
        printf("\nML_SHAKE.c -> MasslessShake() Error: Maximum number of iterations reached!\n");
//        exit(EXIT_FAILURE);
    }
    
    FreeDVector(N, denom);
    FreeDVector(N, sigold);
    
    tfinish = clock();
    
    time = (double) (tfinish - tstart)/CLOCKS_PER_SEC;
    
    PrintStats('S', 0, iter, time, discr, tol);
    WriteStats("logfile.out", "output/", 'S', 0, iter, time, discr, tol);

    return xold;
}
