#ifndef GAUSS_H
#define GAUSS_H

// Метод с двумерным массивом (VLA)
int gauss_2d(int n, double A[n][n], double b[n], double x[n]);

// Метод с одномерным массивом (плоский)
int gauss_1d(int n, double *A, double *b, double *x);

#endif