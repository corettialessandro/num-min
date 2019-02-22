//
//  Constrained.c
//  NumMin
//
//  Created by Alessandro Coretti on 02/22/19.
//  Copyright © 2019 Alessandro Coretti. All rights reserved.
//

#include "Constrained.h"

double * Constrain_b(int N, double ** Ainv, double * b, double x_const, double * constr_b) {

   int i, j;
   double IdAinvb = 0., IdAinvId = 0.;

   for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {

         IdAinvb += Ainv[i][j]*b[j];
         IdAinvId += Ainv[i][j];
      }
   }

   for (i = 0; i < N; i++) {

      constr_b[i] = b[i] + (x_const-IdAinvb)/IdAinvId;
   }

   return constr_b;
}
