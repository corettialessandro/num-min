//
//  ML_BSHAKE.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/21/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "ML_BSHAKE.h"

double * MasslessBlockShake(int N, int nblocks, double ** A, double * b, double * x0, double tol, int maxiter, double * xold) {
    
    clock_t tstart = clock(), tfinish;
    double time;
    
    int iter, block, k, h, kk, hh, i, n, nn;
    double discr = 0;
    double * sigold, * gam;
    double *** Ared, *** Aredinv;
    
    n = (int)nearbyint((double)N/(double)nblocks);
    nn = N%nblocks;
    
    if (nn != 0) {
        
        printf("\nMLB_SHAKE.c -> MasslessBlockShake() Error: Blocks don't divide the matrix exactly. Function not yet implemented!\n");
        exit(EXIT_FAILURE);
    }
    
    sigold = AllocateDVector(N);
    
    Ared = AllocateTensor(nblocks, n, n);
    Aredinv = AllocateTensor(nblocks, n, n);
    gam = AllocateDVector(n);
    
    for (i=0; i<N; i++) xold[i] = x0[i];
    
    for (k=0; k<N; k++) {
        
        sigold[k] = -b[k];
        
        for (i=0; i<N; i++) sigold[k] += A[k][i]*xold[i];
    }
    
    discr = Norm_inf(N, sigold);
    
    for (block=0; block<nblocks; block++) {
        
        hh=0;
        for (h=block*n; h<(block+1)*n; h++) {
            
            kk=0;
            for (k=block*n; k<(block+1)*n; k++) {
                
                Ared[block][hh][kk] = 0.;
                Aredinv[block][hh][kk] = 0.;
                
                for (i=0; i<N; i++) Ared[block][hh][kk] += A[h][i]*A[k][i];
                
                kk++;
            }
            
            hh++;
        }
        
        Aredinv[block] = InvertMatrix(n, Ared[block], Aredinv[block]);
    }
    
    for (iter=0; iter<maxiter; iter++) {
        
        for (block=0; block<nblocks; block++) {
            
            discr = 0;
            
            for (k=block*n; k<(block+1)*n; k++) {
                
                sigold[k] = -b[k];
                
                for (i=0; i<N; i++) sigold[k] += A[k][i]*xold[i];
                
                if (discr < fabs(sigold[k])) discr = fabs(sigold[k]);
            }
            
            if (discr > tol) {
                
                for (kk=0; kk<n; kk++) {
                    
                    gam[kk] = 0;
                    
                    hh=0;
                    for (h=block*n; h<(block+1)*n; h++) {
                        
                        gam[kk] += Aredinv[block][kk][hh]*sigold[h];
                        
                        hh++;
                    }
                }
                
                for (i=0; i<N; i++) {
                    
                    kk=0;
                    for (k=block*n; k<(block+1)*n; k++) {
                        
                        xold[i] -= gam[kk]*A[k][i];
                        
                        kk++;
                    }
                }
            }
        }
        
        discr = Norm_inf(N, sigold);
        
        WriteConvergence('B', "BSH_convergence.out", "output/", iter, discr, maxiter, tol);
        
        if (discr < tol) break;
    }
    
    if (iter == maxiter){
        
        printf("\nMLB_SHAKE.c -> MasslessBlockShake() Error: Maximum number of iterations reached!\n");
//        exit(EXIT_FAILURE);
    }
    
    FreeDVector(N, sigold);
    
    FreeTensor(nblocks, n, n, Ared);
    FreeTensor(nblocks, n, n, Aredinv);
    FreeDVector(n, gam);
    
    tfinish = clock();
    
    time = (double) (tfinish - tstart)/CLOCKS_PER_SEC;
    
    PrintStats('B', nblocks, iter, time, discr, tol);
    WriteStats("logfile.out", "output/", 'B', nblocks, iter, time, discr, tol);
    
    return xold;
}
