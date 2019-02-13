//
//  Tensor.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/20/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "Tensor.h"

double *** AllocateTensor(int rows, int columns, int sections){
    
    double *** tensor;
    
    int row, col;
    
    if ((tensor = (double ***)malloc(rows * sizeof(double **))) == NULL) {
        
        printf("\nTensor.c -> AllocateTensor() Error: Dynamical memory allcocation failed (rows)!\n");
        exit(EXIT_FAILURE);
    }
    
    for (row=0; row<rows; row++) {
        
        if ((tensor[row] = (double **)malloc(columns * sizeof(double *))) == NULL) {
            
            printf("\nTensor.c -> AllocateTensor() Error: Dynamical memory allcocation failed (columns)!\n");
            exit(EXIT_FAILURE);
        }
        
        for (col=0; col<columns; col++) {
                
            if ((tensor[row][col] = (double *)malloc(sections * sizeof(double))) == NULL) {
                    
                printf("\nTensor.c -> AllocateTensor() Error: Dynamical memory allcocation failed (sections)!\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    return tensor;
}

void FreeTensor(int rows, int columns, int sections, double *** tensor) {
    
    int row, col;
    
    for (row = 0; row < rows; row++) {
        
        for (col=0; col<columns; col++) {
            
            free(tensor[row][col]);
        }
        
        free(tensor[row]);
    }
    
    free(tensor);
    
    return;
}
