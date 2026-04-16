#include "gauss.h"
#include <stdio.h>
#include <stdlib.h>

/* Формат входного файла:
 * Строка 1: N
 * Строки 2..N+1: A[i][0] A[i][1] ... A[i][N-1] b[i]
 */
int load_system_from_txt(const char* filename, int* n, double** A, double** b) {
    FILE* f = fopen(filename, "r");
    if (!f) { perror("Ошибка открытия входного файла"); return -1; }

    if (fscanf(f, "%d", n) != 1 || *n <= 0) {
        fprintf(stderr, "Ошибка чтения размера матрицы\n");
        fclose(f); return -1;
    }

    *A = (double*)malloc((*n) * (*n) * sizeof(double));
    *b = (double*)malloc((*n) * sizeof(double));
    if (!*A || !*b) { perror("Ошибка выделения памяти"); free(*A); free(*b); fclose(f); return -1; }

    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
            if (fscanf(f, "%lf", &(*A)[i * (*n) + j]) != 1) goto error;
        }
        if (fscanf(f, "%lf", &(*b)[i]) != 1) goto error;
    }

    fclose(f);
    return 0;

    error:
        free(*A); free(*b); fclose(f);
    fprintf(stderr, "Ошибка чтения данных матрицы\n");
    return -1;
}

int save_solution_to_txt(const char* filename, int n, const double* x) {
    FILE* f = fopen(filename, "w");
    if (!f) { perror("Ошибка открытия файла для записи"); return -1; }

    for (int i = 0; i < n; i++) {
        fprintf(f, "%.10f\n", x[i]);
    }

    fclose(f);
    return 0;
}

// Конвертер 1D -> 2D для совместимости с gauss_2d
double** flat_to_2d(int n, double* A_flat) {
    double** A = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) A[i] = &A_flat[i * n];
    return A;
}