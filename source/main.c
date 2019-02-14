//
//  main.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//
// VER. 0.3

#include <stdio.h>
#include <time.h>

#include "ReadIn.h"
#include "output.h"
#include "Matrix.h"
#include "Vector.h"
#include "ConjGrad.h"
#include "ML_SHAKE.h"
#include "ML_BSHAKE.h"
#include "AdditConstr.h"

#define _MAXSTRLENGTH 50

//TO DO LIST:

int main(int argc, char * argv[]) {

    int N, nblocks, maxiter, dd, verbose = 0, werbose = 0;
    double tol, Akappa, Ar, maxerr_CG, maxerr_SH, maxerr_BSH;
    double * b, * x0, * Adiag, * AR, * Aeigenval, * xF_AN, * xF_CG, * xF_SH, * xF_BSH;
    double ** A, ** Ainv, ** Aeigenvec;
    int additional_constraint;
    double x_const;
    char inputfile[_MAXSTRLENGTH] = "input/input.inpt";
    char Afilename[_MAXSTRLENGTH], bfilename[_MAXSTRLENGTH], x0filename[_MAXSTRLENGTH];
    clock_t an_tstart, an_tfinish;
    double an_time;

    GetOptions(argc, argv, &verbose, &werbose);
    ReadInput(inputfile, &N, &nblocks, Afilename, bfilename, x0filename, &tol, &maxiter, &additional_constraint, &x_const);

    PrintSetup(inputfile, &N, &nblocks, Afilename, bfilename, x0filename, &tol, &maxiter);
    WriteSetup("logfile.out", "output/", &N, &nblocks, Afilename, bfilename, x0filename, &tol, &maxiter);

    A = AllocateMatrix(N, N);
    b = AllocateDVector(N);
    x0 = AllocateDVector(N);
    Adiag = AllocateDVector(N);
    AR = AllocateDVector(N);
    Aeigenvec = AllocateMatrix(N, N);
    Aeigenval = AllocateDVector(N);
    Ainv = AllocateMatrix(N, N);
    xF_AN = AllocateDVector(N);
    xF_CG = AllocateDVector(N);
    xF_SH = AllocateDVector(N);
    xF_BSH = AllocateDVector(N);

    A = ReadMatrix(N, N, A, Afilename);
    b = ReadDVector(N, b, bfilename);
    x0 = ReadDVector(N, x0, x0filename);

    if (verbose) PrintMatrix(N, N, A, "A");
    if (verbose) PrintVector(N, b, "b");
    if (verbose) PrintVector(N, x0, "x0");

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

    if (additional_constraint) {

      ReducedProblem(N, x_const, A, b, x0, verbose, werbose, xF_AN, xF_CG, xF_SH, xF_BSH);

   } else {

      an_tstart = clock();
      Ainv = InvertMatrix(N, A, Ainv);
      xF_AN = RowbyColProd(N, Ainv, b, xF_AN);
      an_tfinish = clock();
      an_time = (double) (an_tfinish - an_tstart)/CLOCKS_PER_SEC;
      PrintStats('A', 0, 0, an_time, 0, 0);
      WriteStats("logfile.out", "output/", 'A', 0, 0, an_time, 0, 0);
      if (verbose) PrintMatrix(N, N, Ainv, "A^{-1}");
      if (verbose) PrintVector(N, xF_AN, "xF_AN");
      if (werbose) WriteMatrix("Ainv.out", "output/", N, N, Ainv, "A^{-1}");
      if (werbose) WriteVector("xF_AN.out", "output/", N, xF_AN, "xF_AN");

      xF_CG = ConjugateGradient(N, A, b, x0, tol, maxiter, xF_CG);
      maxerr_CG = MaxIterError(N, xF_AN, xF_CG);
      WriteError("logfile.out", "output/", 'C', maxerr_CG);
      if (verbose) PrintVector(N, xF_CG, "xF_CG");
      if (werbose) WriteVector("xF_CG.out", "output/", N, xF_CG, "xF_CG");
      printf("  Maximum error component on iterative solution (CG):\n");
      printf("  MAX|xf_AN - xf_CG| = %.4e\n\n", maxerr_CG);

      xF_SH = MasslessShake(N, A, b, x0, tol, maxiter, xF_SH);
      maxerr_SH = MaxIterError(N, xF_AN, xF_SH);
      WriteError("logfile.out", "output/", 'S', maxerr_SH);
      if (verbose) PrintVector(N, xF_SH, "xF_SH");
      if (werbose) WriteVector("xF_SH.out", "output/", N, xF_SH, "xF_SH");
      printf("  Maximum error component on iterative solution (SH):\n ");
      printf("  MAX|xf_AN - xf_SH| = %.4e\n\n", maxerr_SH);

      xF_BSH = MasslessBlockShake(N, nblocks, A, b, x0, tol, maxiter, xF_BSH);
      maxerr_BSH = MaxIterError(N, xF_AN, xF_BSH);
      WriteError("logfile.out", "output/", 'B', maxerr_BSH);
      if (verbose) PrintVector(N, xF_BSH, "xF_BSH");
      if (werbose) WriteVector("xF_BSH.out", "output/", N, xF_BSH, "xF_BSH");
      printf("  Maximum error component on iterative solution (BSH):\n");
      printf("  MAX|xf_AN - xf_BSH| = %.4e\n\n", maxerr_BSH);
   }

    FreeMatrix(N, N, A);
    FreeDVector(N, b);
    FreeDVector(N, x0);
    FreeDVector(N, Adiag);
    FreeDVector(N, AR);
    FreeMatrix(N, N, Aeigenvec);
    FreeDVector(N, Aeigenval);
    FreeMatrix(N, N, Ainv);
    FreeDVector(N, xF_AN);
    FreeDVector(N, xF_CG);
    FreeDVector(N, xF_SH);
    FreeDVector(N, xF_BSH);

    return 0;
}
