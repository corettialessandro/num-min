//
//  Analysis.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/17/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "Analysis.h"

void Analyse(int N, double ** A, double * b, double * x0, int verbose, int werbose) {

   int dd;
   double Akappa, Ar;
   double * AR, * Aeigenval, * Adiag;
   double ** Aeigenvec;

   Adiag = AllocateDVector(N);
   AR = AllocateDVector(N);
   Aeigenvec = AllocateMatrix(N, N);
   Aeigenval = AllocateDVector(N);

   dd = CheckDiagonallyDominance(N, A, Adiag, AR);
   Ar = rcoefficient(N, A);
   if (verbose) PrintVector(N, Adiag, "of diagonal elements of A");
   if (verbose) PrintVector(N, AR, "of sums of off-diagonal elements of A");
   printf("The matrix is");
   dd ? printf(" ") : printf(" not ");
   printf("diagonally dominant!\n");
   printf("r = %.2lf\n\n", Ar);

   Aeigenval = SpectrumMatrix(N, A, Aeigenvec, Aeigenval);
   Akappa = ConditionNumber(N, Aeigenval);
   if (verbose) PrintMatrix(N, N, Aeigenvec, "of eigenvectors of A (column-wise)");
   if (verbose) PrintVector(N, Aeigenval, "of eigenvalues of A");
   printf("Spectral Condition Number (SCN): k = %.4e\n\n", Akappa);
   WriteQFProp("logfile.out", "output/", dd, Ar, Akappa);
   if (werbose) WriteMatrix("Aeigenvec.out", "output/", N, N, Aeigenvec, "of eigenvectors of A (column-wise)");
   if (werbose) WriteVector("Aeigenval.out", "output/", N, Aeigenval, "of eigenvalues of A");

   FreeDVector(N, Adiag);
   FreeDVector(N, AR);
   FreeMatrix(N, N, Aeigenvec);
   FreeDVector(N, Aeigenval);

   return;
}
