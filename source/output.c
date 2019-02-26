//
//  output.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/13/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "output.h"

void PrintSetup(char * inputfile, int dimension, int nconstraints, int nblocks, char * Afilename, char * bfilename, char * x0filename, char * constrfilename, double * constr, double tol, int maxiter) {

   printf("\nNumerical Minimization Test Program\n\n");
   printf("Input Parameters:\n");
   printf("Dimension: N_var = %d\n", dimension);
   printf("Number of constraints: N_constr = %d\n", nconstraints);
   if (nconstraints - dimension == 1) printf("Additional constraint: x_const = %.4lf\n", constr[nconstraints-1]);
   printf("Number of blocks for BSHAKE: nblocks = %d\n", nblocks);
   printf("A Matrix input file: '%s'\n", Afilename);
   printf("B vector input file: '%s'\n", bfilename);
   printf("Initial point input file: '%s'\n", x0filename);
   printf("Values of constraints input file: '%s'\n", constrfilename);
   printf("Tolerance: tol = %.4e\n", tol);
   printf("Maximum number of allowed iterations: maxiter = %d\n", maxiter);
   printf("\n");

   return;
}

void PrintMatrix(int rows, int columns, double ** matrix, char * label) {

    int row, col;

    printf("Matrix %s =\n", label);
    for (row=0; row<rows; row++) {

        printf("|");

        for (col=0; col<columns; col++) {

            printf("%8.4lf\t", matrix[row][col]);
        }

        printf("|\n");
    }

    printf("\n");

    return;
}

void PrintVector(int dim, double * vector, char * label) {

    int row;

    printf("Vector %s =\n", label);

    for (row=0; row<dim; row++) {

        printf("|%8.4lf|\n", vector[row]);
    }

    printf("\n");

    return;
}

void PrintStats(char method, int nblocks, int iter, double exectime, double discr, double tol) {

    int BlockSH = 0, converged = 0, AnSol = 0;
    char meth[100], outcome[100];

    if (method == 'A') {

        sprintf(meth, "ANALYTIC SOLUTION");
        AnSol = 1;

    } else if (method == 'C') {

        sprintf(meth, "CONJUGATE GRADIENT");

    } else if (method == 'S') {

        sprintf(meth, "MASSLESS SHAKE");

    } else if (method == 'B') {

        sprintf(meth, "MASSLESS BLOCK SHAKE");
        BlockSH = 1;

    } else {

        printf("UNRECOGNIZED MINIMIZATION TECHNIQUE!\n");
        exit(EXIT_FAILURE);
    }

    if (discr <= tol) {

        sprintf(outcome, "SUCCEEDED");
        converged = 1;

    } else {

        sprintf(outcome, "FAILED");
    }

    printf("\n");

    printf("* STATISTICS FOR %s (%s):\n", meth, outcome);
    if (BlockSH) printf("  NUMBER OF BLOCKS = %d\n", nblocks);
    printf("  TIME OF EXECUTION = %.4e s\n", exectime);
    if (!AnSol) {
        if (converged) {
            printf("  NUMBER OF ITERATIONS FOR CONVERGENCE = %d\n", iter);
            printf("  DISCRIMINANT AT CONVERGENCE = %.4e\n", discr);
        } else {
            printf("  MAXIMUM NUMBER OF ITERATIONS REACHED: maxiter = %d\n", iter);
            printf("  DISCRIMINANT AT 'maxiter' = %.4e\n", discr);
        }
        printf("  TOLERANCE = %.4e\n", tol);
    }
    printf("\n");

    return;
}

void WriteSetup(char * outputfile, char * outputpath, int dimension, int nconstraints, int nblocks, char * Afilename, char * bfilename, char * x0filename, char * constrfilename, double * constr, double tol, int maxiter) {

    FILE * fp_setup;

    char output[100];

    sprintf(output, "%s%s", outputpath, outputfile);

    if ((fp_setup = fopen(output, "w+")) == NULL){
        printf("\noutput.c -> WriteSetup() Error: File '%s' not found!\n", output);
        exit(EXIT_FAILURE);
   }

   fprintf(fp_setup, "\nNumerical Minimization Test Program\n\n");
   fprintf(fp_setup, "Input Parameters:\n");
   fprintf(fp_setup, "Dimension: N_var = %d\n", dimension);
   fprintf(fp_setup, "Number of constraints: N_constr = %d\n", nconstraints);
   if (nconstraints - dimension == 1) fprintf(fp_setup, "Additional constraint: x_const = %.4lf\n", constr[nconstraints-1]);
   fprintf(fp_setup, "Number of blocks for BSHAKE: nblocks = %d\n", nblocks);
   fprintf(fp_setup, "A Matrix input file: '%s'\n", Afilename);
   fprintf(fp_setup, "B vector input file: '%s'\n", bfilename);
   fprintf(fp_setup, "Initial point input file: '%s'\n", x0filename);
   fprintf(fp_setup, "Values of constraints input file: '%s'\n", constrfilename);
   fprintf(fp_setup, "Tolerance: tol = %.4e\n", tol);
   fprintf(fp_setup, "Maximum number of allowed iterations: maxiter = %d\n", maxiter);

   fclose(fp_setup);

   return;
}

void WriteMatrix(char * outputfile, char * outputpath, int rows, int columns, double ** matrix, char * label) {

    FILE * fp_matrix;

    char output[100];

    sprintf(output, "%s%s", outputpath, outputfile);

    if ((fp_matrix = fopen(output, "w+")) == NULL){
        printf("\noutput.c -> WriteMatrix() Error: File '%s' not found!\n", output);
        exit(EXIT_FAILURE);
    }

    int row, col;

    fprintf(fp_matrix, "# Matrix %s =\n", label);
    for (row=0; row<rows; row++) {

        fprintf(fp_matrix, "|");

        for (col=0; col<columns; col++) {

            fprintf(fp_matrix, "%8.4e\t", matrix[row][col]);
        }

        fprintf(fp_matrix, "|\n");
    }

    fprintf(fp_matrix, "\n");

    fclose(fp_matrix);

    return;
}

void WriteVector(char * outputfile, char * outputpath, int dim, double * vector, char * label) {

    FILE * fp_vector;

    char output[100];

    sprintf(output, "%s%s", outputpath, outputfile);

    if ((fp_vector = fopen(output, "w+")) == NULL){
        printf("\noutput.c -> WriteMatrix() Error: File '%s' not found!\n", output);
        exit(EXIT_FAILURE);
    }

    int row;

    fprintf(fp_vector, "# Vector %s =\n", label);

    for (row=0; row<dim; row++) {

        fprintf(fp_vector, "|%8.4e|\n", vector[row]);
    }

    fprintf(fp_vector, "\n");

    return;
}

void WriteQFProp(char * outputfile, char * outputpath, int dd, double r, double scn){

    FILE *fp_QFProp;

    char output[100];

    sprintf(output, "%s%s", outputpath, outputfile);

    if ((fp_QFProp = fopen(output, "a")) == NULL){
        printf("\noutput.c -> WriteStats() Error: File '%s' not found!\n", output);
        exit(EXIT_FAILURE);
    }

    fprintf(fp_QFProp, "\nThe matrix is");
    dd ? fprintf(fp_QFProp, " ") : fprintf(fp_QFProp, " not ");
    fprintf(fp_QFProp, "diagonally dominant!\n");
    fprintf(fp_QFProp, "r = %.2lf\n\n", r);
    fprintf(fp_QFProp, "Spectral Condition Number (SCN): k = %.4e\n\n", scn);

    fclose(fp_QFProp);

    return;
}

void WriteStats(char * outputfile, char * outputpath, char method, int nblocks, int iter, double exectime, double discr, double tol) {

    FILE *fp_stats;

    char output[100];

    sprintf(output, "%s%s", outputpath, outputfile);

    if ((fp_stats = fopen(output, "a")) == NULL){
        printf("\noutput.c -> WriteStats() Error: File '%s' not found!\n", output);
        exit(EXIT_FAILURE);
    }

    int BlockSH = 0, converged = 0, AnSol = 0;
    char meth[100], outcome[100];

    if (method == 'A') {

        sprintf(meth, "ANALYTIC SOLUTION");
        AnSol = 1;

    } else if (method == 'C') {

        sprintf(meth, "CONJUGATE GRADIENT");

    } else if (method == 'S') {

        sprintf(meth, "MASSLESS SHAKE");

    } else if (method == 'B') {

        sprintf(meth, "MASSLESS BLOCK SHAKE");
        BlockSH = 1;

    } else {

        printf("UNRECOGNIZED MINIMIZATION TECHNIQUE!\n");
        exit(EXIT_FAILURE);
    }

    if (discr <= tol) {

        sprintf(outcome, "SUCCEEDED");
        converged = 1;

    } else {

        sprintf(outcome, "FAILED");
    }

    fprintf(fp_stats, "\n");

    fprintf(fp_stats, "* STATISTICS FOR %s (%s):\n", meth, outcome);
    if (BlockSH) fprintf(fp_stats, "  NUMBER OF BLOCKS = %d\n", nblocks);
    fprintf(fp_stats, "  TIME OF EXECUTION = %.4e s\n", exectime);
    if (!AnSol) {
        if (converged) {
            fprintf(fp_stats, "  NUMBER OF ITERATIONS FOR CONVERGENCE = %d\n", iter);
            fprintf(fp_stats, "  DISCRIMINANT AT CONVERGENCE = %.4e\n", discr);
        } else {
            fprintf(fp_stats, "  MAXIMUM NUMBER OF ITERATIONS REACHED: maxiter = %d\n", iter);
            fprintf(fp_stats, "  DISCRIMINANT AT 'maxiter' = %.4e\n", discr);
        }
        fprintf(fp_stats, "  TOLERANCE = %.4e\n", tol);
    }
    fprintf(fp_stats, "\n");

    fclose(fp_stats);

    return;
}

void WriteError(char * outputfile, char * outputpath, char method, double maxerr) {

    FILE *fp_error;

    char output[100], meth[100];

    sprintf(output, "%s%s", outputpath, outputfile);

    if ((fp_error = fopen(output, "a")) == NULL){
        printf("\noutput.c -> WriteError() Error: File '%s' not found!\n", output);
        exit(EXIT_FAILURE);
    }

    if (method == 'C') {

        sprintf(meth, "(CG):\n  MAX|xf_AN - xf_CG|");

    } else if (method == 'S') {

        sprintf(meth, "(SH):\n  MAX|xf_AN - xf_SH|");

    } else if (method == 'B') {

        sprintf(meth, "(BSH):\n  MAX|xf_AN - xf_BSH|");

    } else {

        printf("UNRECOGNIZED MINIMIZATION TECHNIQUE!\n");
        exit(EXIT_FAILURE);
    }

    fprintf(fp_error, "  Maximum error component on iterative solution %s = %.4e\n\n", meth, maxerr);

    fclose(fp_error);

    return;
}

void WriteConvergence(char method, char * outputfile, char * outputpath, int iter, double discr, int maxiter, double tol) {

    FILE * fp_convergence;

    char output[100];

    sprintf(output, "%s%s", outputpath, outputfile);

    if (iter == 0) {

        if ((fp_convergence = fopen(output, "w+")) == NULL){
            printf("\noutput.c -> WriteConvergence() Error: File '%s' not found!\n", output);
            exit(EXIT_FAILURE);
        }

        if (method == 'C') {

            fprintf(fp_convergence, "# MINIMIZATION PERFORMED WITH CONJUGATE GRADIENT METHOD\n");

        } else if (method == 'S') {

            fprintf(fp_convergence, "# MINIMIZATION PERFORMED WITH SHAKE METHOD\n");
        }

        fprintf(fp_convergence, "# MAXITER = %d\n", maxiter);
        fprintf(fp_convergence, "# TOL = %.4e\n", tol);

        fprintf(fp_convergence, "# iterations\tdiscriminant\n");
        fprintf(fp_convergence, "%d\t%.4e\n", iter, discr);

    } else {

        if ((fp_convergence = fopen(output, "a")) == NULL){
            printf("\noutput.c -> WriteConvergence() Error: File '%s' not found!\n", output);
            exit(EXIT_FAILURE);
        }

        fprintf(fp_convergence, "%d\t%.4e\n", iter, discr);
    }

    fclose(fp_convergence);

    return;
}
