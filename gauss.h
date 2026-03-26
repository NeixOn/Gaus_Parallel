#ifndef GAUSS_H
#define GAUSS_H
#include <stdbool.h>

// Последовательные версии - используем указатели вместо VLA
double* gauss_2d(int n, double **A, double *b);   // 2D: массив указателей
double* gauss_1d(int n, double *A, double *b);    // 1D: плоский массив

// Параллельные версии
int gauss_2d_parallel(int N);
int gauss_1d_parallel(int N);

// Тесты
void run_test(int n, bool show_matrices);

#endif