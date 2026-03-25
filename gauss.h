#ifndef GAUSS_H
#define GAUSS_H

// Метод с двумерным массивом (VLA)
double* gauss_2d(int n, double A[n][n], double b[n]);

// Метод с одномерным массивом (плоский)
double* gauss_1d(int n, double *A, double *b);

#endif