#ifndef GAUSS_H
#define GAUSS_H
#include <stdbool.h>

// Метод с двумерным массивом (VLA)
double* gauss_2d(int n, double A[n][n], double b[n]);

int gauss_2d_parallel(int N);

int gauss_1d_parallel(int N);

// Метод с одномерным массивом (плоский)
double* gauss_1d(int n, double *A, double *b);

void run_test(int n, bool show_matrices);

#endif