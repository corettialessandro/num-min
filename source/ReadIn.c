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
            case 'v':
                * verbose = 1;
                break;
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

void ReadInput(char * inputfile, int * dimension, int * nblocks, char * Afilename, char * bfilename, char * x0filename, double * tol, int * maxiter) {
    
    FILE *fp_input;
    
    if ((fp_input = fopen(inputfile, "r+")) == NULL){
        printf("\nReadIn.c -> ReadInput() Error: File '%s' not found!\n", inputfile);
        exit(EXIT_FAILURE);
    }
    
    fscanf(fp_input, "%d %*[^\n]", dimension);
    fscanf(fp_input, "%d %*[^\n]", nblocks);
    fscanf(fp_input, "%s %*[^\n]", Afilename);
    fscanf(fp_input, "%s %*[^\n]", bfilename);
    fscanf(fp_input, "%s %*[^\n]", x0filename);
    fscanf(fp_input, "%lf %*[^\n]", tol);
    fscanf(fp_input, "%d %*[^\n]", maxiter);
    
    fclose(fp_input);
    
    return;
}


