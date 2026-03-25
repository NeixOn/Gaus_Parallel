#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gauss.h"

#define N 2500  // Размер системы

// Выбираем метод: 1 - двумерный, 2 - одномерный
#define METHOD 1

int main() {
    system("chcp 65001");
    // Для двумерного метода
#if METHOD == 1
    double (*A)[N] = malloc(sizeof(double[N][N]));
    double *b = malloc(N * sizeof(double));
    double *x = malloc(N * sizeof(double));

    if (!A || !b || !x) {
        printf("Ошибка выделения памяти\n");
        return 1;
    }

    // Заполняем матрицу (диагональное преобладание)
    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        double rowSum = 0;
        for (int j = 0; j < N; j++) {
            if (i == j) {
                A[i][j] = (double)(rand() % 100 + 50);
            } else {
                A[i][j] = (double)(rand() % 20 - 10);
                rowSum += fabs(A[i][j]);
            }
        }
        A[i][i] += rowSum + 1;
        b[i] = (double)(rand() % 200 - 100);
    }

    printf("Система %d x %d создана (двумерный массив)\n", N, N);

    clock_t start = clock();
    int success = gauss_2d(N, A, b, x);
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

#endif

    // Для одномерного метода
#if METHOD == 2
    double *A = malloc(N * N * sizeof(double));
    double *b = malloc(N * sizeof(double));
    double *x = malloc(N * sizeof(double));

    if (!A || !b || !x) {
        printf("Ошибка выделения памяти\n");
        return 1;
    }

    // Заполняем матрицу (диагональное преобладание)
    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        double rowSum = 0;
        for (int j = 0; j < N; j++) {
            if (i == j) {
                A[i * N + j] = (double)(rand() % 100 + 50);
            } else {
                A[i * N + j] = (double)(rand() % 20 - 10);
                rowSum += fabs(A[i * N + j]);
            }
        }
        A[i * N + i] += rowSum + 1;
        b[i] = (double)(rand() % 200 - 100);
    }

    printf("Система %d x %d создана (одномерный массив)\n", N, N);

    clock_t start = clock();
    int success = gauss_1d(N, A, b, x);
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

#endif

    if (success) {
        printf("Метод Гаусса выполнен успешно\n");
        printf("Время выполнения: %.4f секунд\n", time_spent);

        // Выводим первые 5 решений
        printf("\nПервые 5 решений:\n");
        for (int i = 0; i < (N < 5 ? N : 5); i++) {
            printf("x[%d] = %.6f\n", i, x[i]);
        }

        // Проверка
        printf("\nПроверка (Ax - b) для первых 5 уравнений:\n");
#if METHOD == 1
        for (int i = 0; i < (N < 5 ? N : 5); i++) {
            double sum = 0;
            for (int j = 0; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            printf("|%d: %.2e|\n", i, fabs(sum - b[i]));
        }
#endif

#if METHOD == 2
        for (int i = 0; i < (N < 5 ? N : 5); i++) {
            double sum = 0;
            for (int j = 0; j < N; j++) {
                sum += A[i * N + j] * x[j];
            }
            printf("|%d: %.2e|\n", i, fabs(sum - b[i]));
        }
#endif
    } else {
        printf("Ошибка при решении системы\n");
    }

    free(A);
    free(b);
    free(x);

    return 0;
}