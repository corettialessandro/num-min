//
//  Matrix.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "Matrix.h"

double ** AllocateMatrix(int rows, int columns) {
    
    double ** matrix;
    
    if ((matrix = (double **)malloc(rows * sizeof(double *))) == NULL) {
        
        printf("\nMatrix.c -> AllocateMatrix() Error: Dynamical memory allocation failed (rows)!\n");
        exit(EXIT_FAILURE);
    }
    
    int row;
    
    for (row=0; row<rows; row++) {
        
        if ((matrix[row] = (double *)malloc(columns * sizeof(double))) == NULL) {

            printf("\nMatrix.c -> AllocateMatrix() Error: Dynamical memory allocation failed (columns)!\n");
            exit(EXIT_FAILURE);
        }
    }
    
    return matrix;
}

double ** ReadMatrix(int rows, int columns, double ** matrix, char * filepath) {
    
    FILE * fp_matrix;
    
    int row, col = 0;
    
    if ((fp_matrix = fopen(filepath, "r+")) == NULL){
        printf("\nMatrix.c -> ReadMatrix() Error: File '%s' not found!\n", filepath);
        exit(EXIT_FAILURE);
    }

    for (row=0; row<rows; row++) {
        
        for (col=0; col<columns; col++) {
            
            if ((fscanf(fp_matrix, "%lf", &matrix[row][col])) <= 0) {
             
                printf("\nMatrix.c -> ReadMatrix() Error on (%d, %d): Wrong matrix dimension!\n", row, col);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    fclose(fp_matrix);
    
    return matrix;
}

double ** InvertMatrix(int dimension, double ** matrix, double ** inverse) {
    
    int i, j, info, lwork = dimension*dimension;
    int * ipiv;
    double * work, * mattmp;

    ipiv = AllocateIVector(dimension);
    work = AllocateDVector(dimension*dimension);
    mattmp = AllocateDVector(dimension*dimension);
    
    for (i=0; i<dimension; i++) {
        
        for (j=0; j<dimension; j++) {
            
            mattmp[i*dimension + j] = matrix[i][j];
        }
    }
    
    dgetrf_(&dimension, &dimension, mattmp, &dimension, ipiv, &info);

    if (info != 0) {
        
        printf("\nMatrix.c -> InvertMatrix() Error: Matrix is numerically singular!\n");
        exit(EXIT_FAILURE);
    }
    
    dgetri_(&dimension, mattmp, &dimension, ipiv, work, &lwork, &info);

    if (info != 0) {
        
        printf("\nMatrix.c -> InvertMatrix() Error: Matrix inversion failed!\n");
        exit(EXIT_FAILURE);
    }

    for (i=0; i<dimension; i++) {
        
        for (j=0; j<dimension; j++) {
            
            inverse[i][j] = mattmp[i*dimension + j];
        }
    }

    FreeIVector(dimension, ipiv);
    FreeDVector(dimension*dimension, work);
    FreeDVector(dimension*dimension, mattmp);
    
    return inverse;
}

double * SpectrumMatrix(int dimension, double ** matrix, double ** eigen_matrix, double * eigenvalues){
    
    int i, j;
    int info, lwork = -1, lda = dimension;
    double wkopt;
    double * work, * mattmp;
    
    mattmp = AllocateDVector(dimension*dimension);
    
    for (i=0; i<dimension; i++) {
        
        for (j=0; j<dimension; j++) {
            
            mattmp[i*dimension + j] = matrix[i][j];
        }
    }

    dsyev_("Vectors", "Upper", &dimension, mattmp, &lda, eigenvalues, &wkopt, &lwork, &info);

    lwork = (int)wkopt;

    work = AllocateDVector(lwork);
    
    dsyev_("Vectors", "Upper", &dimension, mattmp, &lda, eigenvalues, work, &lwork, &info);
    
    if (info != 0) {
        
        printf("\nMatrix.c -> InvertMatrix() Error: Matrix eigenvalue computation failed!\n");
        exit(EXIT_FAILURE);
    }

    for (i=0; i<dimension; i++) {
        
        for (j=0; j<dimension; j++) {
            
            eigen_matrix[i][j] = mattmp[i*dimension + j];
        }
    }

    FreeDVector(lwork, work);
    FreeDVector(dimension*dimension, mattmp);

    return eigenvalues;
}

double ConditionNumber(int dimension, double * eigenvalues) {
    
    double kappa;
    
    int i;
    double lmin = eigenvalues[0], lmax = eigenvalues[dimension-1];
    
    for (i=0; i<dimension; i++) {
        
        if (lmin > eigenvalues[i]) {
            
            lmin = eigenvalues[i];
        }
        
        if (lmax < eigenvalues[i]) {
            
            lmax = eigenvalues[i];
        }
    }
    
    kappa = lmax/lmin;
    
    return kappa;
}

int CheckDiagonallyDominance(int dimension, double ** matrix, double * diagonal, double * R){
    
    int i, j, counter = 0;
    
    for (i=0; i<dimension; i++) {
        
        R[i] = 0;
        
        for (j=0; j<dimension; j++) {
            
            if (i!=j) R[i] += fabs(matrix[i][j]);
        }
        
        diagonal[i] = fabs(matrix[i][i]);
        
        if (diagonal[i] >= R[i]) counter++;
    }
    
    if (counter == dimension) {
        
        return 1;
    
    } else {
        
        return 0;
    }
}

double rcoefficient(int dimension, double ** matrix) {
    
    double r;
    
    int i, j, ip1, jp1;
    double absmij, n = 0, sumx = 0, sumy = 0, sumxsq = 0, sumysq = 0, sumxy = 0;
    
    for (i=0; i<dimension; i++) {
        
        ip1 = i+1;
        
        for (j=0; j<dimension; j++) {
        
            jp1 = j+1;
            
            absmij = fabs(matrix[i][j]);
            
            n += absmij;
            sumx += ip1*absmij;
            sumy += jp1*absmij;
            sumxsq += ip1*ip1*absmij;
            sumysq += jp1*jp1*absmij;
            sumxy += ip1*jp1*absmij;
        }
    }
    
    r = (n*sumxy - sumx*sumy)/(sqrt(n*sumxsq - sumx*sumx)*sqrt(n*sumysq - sumy*sumy));

    return r;
}

void FreeMatrix(int rows, int columns, double ** matrix) {
    
    int row;
    
    for (row = 0; row < rows; row++) {
        
        free(matrix[row]);
    }
    
    free(matrix);
    
    return;
}
