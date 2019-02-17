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
#include "Analysis.h"
#include "Minimization.h"

#define _MAXSTRLENGTH 50

//TO DO LIST:

int main(int argc, char * argv[]) {

    int N, nblocks, maxiter, verbose = 0, werbose = 0;
    double tol;
    double * b, * x0, * xF_AN, * xF_CG, * xF_SH, * xF_BSH;
    double ** A;
    int constrained;
    double x_const;
    char inputfile[_MAXSTRLENGTH] = "input/input.inpt";
    char Afilename[_MAXSTRLENGTH], bfilename[_MAXSTRLENGTH], x0filename[_MAXSTRLENGTH];

    GetOptions(argc, argv, &verbose, &werbose);
    ReadInput(inputfile, &N, &nblocks, Afilename, bfilename, x0filename, &tol, &maxiter, &constrained, &x_const);

    PrintSetup(inputfile, N, nblocks, Afilename, bfilename, x0filename, tol, maxiter, constrained, x_const);
    WriteSetup("logfile.out", "output/", N, nblocks, Afilename, bfilename, x0filename, tol, maxiter, constrained, x_const);

    A = AllocateMatrix(N, N);
    b = AllocateDVector(N);
    x0 = AllocateDVector(N);
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

    Analyse(N, A, b, x0, verbose, werbose);

    if (constrained) {

      Reduced(N, x_const, A, b, x0, verbose, werbose, tol, maxiter, nblocks, xF_AN, xF_CG, xF_SH, xF_BSH);

   } else {

      Unconstrained(N, A, b, x0, verbose, werbose, tol, maxiter, nblocks, xF_AN, xF_CG, xF_SH, xF_BSH);
   }

    FreeMatrix(N, N, A);
    FreeDVector(N, b);
    FreeDVector(N, x0);
    FreeDVector(N, xF_AN);
    FreeDVector(N, xF_CG);
    FreeDVector(N, xF_SH);
    FreeDVector(N, xF_BSH);

    return 0;
}
