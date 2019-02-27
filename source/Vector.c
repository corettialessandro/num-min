//
//  Vector.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "Vector.h"

double * AllocateDVector(int dim) {

    double * vector;

    if ((vector = (double *)malloc(dim * sizeof(double))) == NULL) {

        printf("\nVector.c -> AllocateDVector() Error: Dynamical memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    return vector;
}

int * AllocateIVector(int dim) {

    int * vector;

    if ((vector = (int *)malloc(dim * sizeof(int))) == NULL) {

        printf("\nVector.c -> AllocateIVector() Error: Dynamical memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    return vector;
}

double * CAllocateDVector(int dim) {

    double * vector;

    if ((vector = (double *)calloc(dim, sizeof(double))) == NULL) {

        printf("\nVector.c -> CAllocateDVector() Error: Dynamical memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    return vector;
}

double * ReadDVector(int dim, double * vector, char * filepath) {

    FILE * fp_vector;

    int row;

    if ((fp_vector = fopen(filepath, "r+")) == NULL){
        printf("\nVector.c -> ReadDVector() Error: File '%s' not found!\n", filepath);
        exit(EXIT_FAILURE);
    }

    for (row=0; row<dim; row++) {

        if (fscanf(fp_vector, "%lf", &vector[row]) <= 0) {

            printf("\nVector.c -> ReadDVector() Error on element %d: Wrong vector dimension!\n", row);
            exit(EXIT_FAILURE);
        }
    }

    fclose(fp_vector);

    return vector;
}

void FreeDVector(int dim, double * vector) {

    free(vector);

    return;
}

void FreeIVector(int dim, int * vector) {

    free(vector);

    return;
}
