#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include "gauss.h"

// Создаёт случайную матрицу с диагональным преобладанием
// Теперь принимает double ** вместо VLA
void create_system(int n, double **A2D, double *A1D, double *b) {
    for (int i = 0; i < n; i++) {
        double rowSum = 0;
        for (int j = 0; j < n; j++) {
            double val;
            if (i == j) {
                val = (double)(rand() % 100 + 50);
            } else {
                val = (double)(rand() % 20 - 10);
                rowSum += fabs(val);
            }
            A2D[i][j] = val;
            A1D[i * n + j] = val;
        }
        A2D[i][i] += rowSum + 1;
        A1D[i * n + i] += rowSum + 1;
        b[i] = (double)(rand() % 200 - 100);
    }
}

// Тестирует метод Гаусса для заданного N
void run_test(int n, bool show_matrices) {
    // Выделяем память: 2D массив как массив указателей
    double **A2D = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A2D[i] = (double *)malloc(n * sizeof(double));
    }
    double *A1D = malloc(n * n * sizeof(double));
    double *b = malloc(n * sizeof(double));
    double *x2 = NULL;
    double *x1 = NULL;

    if (!A2D || !A1D || !b) {
        free(A1D);
        free(b);
        for (int i = 0; i < n; i++) free(A2D[i]);
        free(A2D);
        return;
    }

    srand(time(NULL));
    create_system(n, A2D, A1D, b);

    // Замер для 2D метода
    clock_t start = clock();
    x2 = gauss_2d(n, A2D, b);
    clock_t end = clock();
    double time_2d = (double)(end - start) / CLOCKS_PER_SEC;

    // Замер для 1D метода
    start = clock();
    x1 = gauss_1d(n, A1D, b);
    end = clock();
    double time_1d = (double)(end - start) / CLOCKS_PER_SEC;

    // Освобождаем память
    for (int i = 0; i < n; i++) free(A2D[i]);
    free(A2D);
    free(A1D);
    free(b);
    if (x2) free(x2);
    if (x1) free(x1);

    // Вывод результатов
    printf("\n========================================\n");
    printf("size matrix: %d\n", n);
    printf("Results:\n");
    printf("2d: %.4f second\n", time_2d);
    printf("1d: %.4f second\n", time_1d);
}