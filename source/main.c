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
#include "Minimization.h"

#define _MAXSTRLENGTH 50

//TO DO LIST:

int main(int argc, char * argv[]) {

    int N_var, N_constr, nblocks, maxiter, verbose = 0, werbose = 0;
    double tol;
    double * b, * x0, * constr, * xF_AN, * xF_CG, * xF_SH, * xF_BSH;
    double ** A;
    char inputfile[_MAXSTRLENGTH] = "input/input.inpt";
    char Afilename[_MAXSTRLENGTH], bfilename[_MAXSTRLENGTH], x0filename[_MAXSTRLENGTH], constrfilename[_MAXSTRLENGTH];

    GetOptions(argc, argv, &verbose, &werbose);
    ReadInput(inputfile, &N_var, &N_constr, &nblocks, Afilename, bfilename, x0filename, constrfilename, &tol, &maxiter);

    PrintSetup(inputfile, N_var, N_constr, nblocks, Afilename, bfilename, x0filename, constrfilename, tol, maxiter, 0);
    WriteSetup("logfile.out", "output/", N_var, N_constr, nblocks, Afilename, bfilename, x0filename, constrfilename, tol, maxiter, 0);

    Initialize(N_var, N_constr, &A, &b, &x0, &constr, Afilename, bfilename, x0filename, constrfilename, verbose, &xF_AN, &xF_CG, &xF_SH, &xF_BSH);

    Analyse(N_var, A, b, x0, verbose, werbose);

    Minimize(N_var, N_constr, A, b, x0, constr, verbose, werbose, tol, maxiter, nblocks, xF_AN, xF_CG, xF_SH, xF_BSH);

    Finalize(N_var, N_constr, &A, &b, &x0, &constr, &xF_AN, &xF_CG, &xF_SH, &xF_BSH);

    return 0;
}
