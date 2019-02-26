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
#include "setup.h"
#include "output.h"
#include "Matrix.h"
#include "Vector.h"
#include "Analysis.h"
#include "Unconstrained.h"
#include "Reduced.h"
#include "Constrained.h"

#define _MAXSTRLENGTH 50

//TO DO LIST:

int main(int argc, char * argv[]) {

    int N, nblocks, maxiter, verbose = 0, werbose = 0;
    double tol;
    double * b, * x0, * constr, * xF_AN, * xF_CG, * xF_SH, * xF_BSH;
    double ** A;
    char mode;
    double x_const = 0;
    char inputfile[_MAXSTRLENGTH] = "input/input.inpt";
    char Afilename[_MAXSTRLENGTH], bfilename[_MAXSTRLENGTH], x0filename[_MAXSTRLENGTH], constrfilename[_MAXSTRLENGTH];

    GetOptions(argc, argv, &verbose, &werbose);
    ReadInput(inputfile, &N, &nblocks, Afilename, bfilename, x0filename, constrfilename, &tol, &maxiter, &mode, &x_const);

    PrintSetup(inputfile, N, nblocks, Afilename, bfilename, x0filename, constrfilename, tol, maxiter, mode, x_const);
    WriteSetup("logfile.out", "output/", N, nblocks, Afilename, bfilename, x0filename, constrfilename, tol, maxiter, mode, x_const);

    Initialize(N, &A, &b, &x0, &constr, Afilename, bfilename, x0filename, constrfilename, verbose, &xF_AN, &xF_CG, &xF_SH, &xF_BSH);

    Analyse(N, A, b, x0, verbose, werbose);

   switch (mode) {

      case 'U':
         Unconstrained(N, A, b, x0, constr, verbose, werbose, tol, maxiter, nblocks, xF_AN, xF_CG, xF_SH, xF_BSH);
         break;
      case 'R':
         Reduced(N, x_const, A, b, x0, constr, verbose, werbose, tol, maxiter, nblocks, xF_AN, xF_CG, xF_SH, xF_BSH);
         break;
      case 'C':
         Constrained(N, A, b, x0, constr, verbose, werbose, tol, maxiter, x_const, xF_AN, xF_SH);
         break;
      default:
         printf("Unrecognized execution mode: mode = %c.\nTerminating.\n\n", mode);
         exit(EXIT_FAILURE);
   }

    Finalize(N, &A, &b, &x0, &constr, &xF_AN, &xF_CG, &xF_SH, &xF_BSH);

    return 0;
}
