//
//  aux_func.c
//  NumMin
//
//  Created by Alessandro Coretti on 12/12/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "aux_func.h"

double Norm_inf(int N, double * v) {

    double norm = 0;

    int i;

    for (i=0; i<N; i++) {

        if (fabs(v[i]) > norm) norm = fabs(v[i]);
    }

    return norm;
}

double Norm_l2(int N, double * v) {

    double norm = 0;

    int i;

    for (i=0; i<N; i++) {

        norm += v[i]*v[i];
    }

    norm = sqrt(norm);

    return norm;
}

double modsq(int N, double * v) {

    double modsq = 0;

    int i;

    for (i=0; i<N; i++) {

        modsq += v[i]*v[i];
    }

    return modsq;
}

double ScalProd(int N, double * u, double * v) {

    double udotv = 0;

    int i;

    for (i=0; i<N; i++) {

        udotv += u[i]*v[i];
    }

    return udotv;
}

double * RowbyColProd(int N, double ** M, double * v, double * result) {

    int i, j;

    for (i=0; i<N; i++) {

        result[i] = 0;

        for (j=0; j<N; j++) {

            result[i] += M[i][j]*v[j];
        }
    }

    return result;
}

double * VectorSum(int N, double * u, double * v, double * result) {

    int i;

    for (i=0; i<N; i++) {

        result[i] = u[i] + v[i];
    }

    return result;
}

double * VectorDifference(int N, double * u, double * v, double * result) {

    int i;

    for (i=0; i<N; i++) {

        result[i] = u[i] - v[i];
    }

    return result;
}

double * VectorScalarSum(int N, double * u, double l, double * result) {

   int i;

   for (i = 0; i < N; i++) {

      result[i] = u[i] + l;
   }

   return result;
}

double MaxIterError(int N, double * an_sol, double * iter_sol) {

    double maxerr;

    int i;
    double * error;

    error = AllocateDVector(N);

    for (i=0; i<N; i++) {

        error[i] = an_sol[i] - iter_sol[i];
    }

    maxerr = Norm_inf(N, error);

    FreeDVector(N, error);

    return maxerr;
}

double SumComponents(int N, double * x) {

   double sum_components = 0.;

   int i;

   for (i = 0; i < N; i++) {

      sum_components += x[i];
   }

   return sum_components;
}
