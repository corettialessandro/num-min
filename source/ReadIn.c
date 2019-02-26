//
//  ReadIn.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/13/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "ReadIn.h"

void GetOptions(int argc, char * argv[], int * verbose, int * werbose) {

    int c;

    while ((c = getopt(argc, argv, "vw")) != -1) {
        switch (c) {
            // verbose on screen
            case 'v':
                * verbose = 1;
                break;
            // verbose on files
            case 'w':
                * werbose = 1;
                break;
            case '?':
                if (isprint(optopt)) {
                    printf("Unknown option '-%c'.\nExecution aborted.\n\n", optopt);
                } else {
                    printf("Unknown option character '\\x%x'.Execution aborted.\n\n", optopt);
                }
                exit(EXIT_FAILURE);
            default:
                abort();
        }
    }

    return;
}

void ReadInput(char * inputfile, int * dimension, int * nconstraints, int * nblocks, char * Afilename, char * bfilename, char * x0filename, char * constrfilename, double * tol, int * maxiter, char * mode, double * x_const) {

    FILE *fp_input;

    if ((fp_input = fopen(inputfile, "r+")) == NULL){
        printf("\nReadIn.c -> ReadInput() Error: File '%s' not found!\n", inputfile);
        exit(EXIT_FAILURE);
    }

    fscanf(fp_input, "%d %*[^\n]", dimension);
    fscanf(fp_input, "%d %*[^\n]", nconstraints);
    fscanf(fp_input, "%d %*[^\n]", nblocks);
    fscanf(fp_input, "%s %*[^\n]", Afilename);
    fscanf(fp_input, "%s %*[^\n]", bfilename);
    fscanf(fp_input, "%s %*[^\n]", x0filename);
    fscanf(fp_input, "%s %*[^\n]", constrfilename);
    fscanf(fp_input, "%lf %*[^\n]", tol);
    fscanf(fp_input, "%d %*[^\n]", maxiter);
    fscanf(fp_input, " %c %*[^\n]", mode);
    if (*mode != 'U') fscanf(fp_input, "%lf %*[^\n]", x_const);

    fclose(fp_input);

    return;
}
