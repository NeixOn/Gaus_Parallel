#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include "gauss.h"

#define N 2000  // Для теста маленький размер, чтобы видеть вывод

// Флаг вывода матриц (true - выводить, false - не выводить)
#define SHOW_MATRICES false

int main() {
    system("chcp 65001");

    // Выделяем память
    double (*A2D)[N] = malloc(sizeof(double[N][N]));
    double *A1D = malloc(N * N * sizeof(double));
    double *b = malloc(N * sizeof(double));
    double *x2 = NULL;
    double *x1 = NULL;

    if (!A2D || !A1D || !b) {
        printf("Ошибка выделения памяти\n");
        free(A2D);
        free(A1D);
        free(b);
        return 1;
    }

    // СОЗДАЁМ ОДНУ МАТРИЦУ
    printf("Создание системы %d x %d...\n", N, N);
    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        double rowSum = 0;
        for (int j = 0; j < N; j++) {
            double val;
            if (i == j) {
                val = (double)(rand() % 100 + 50);
            } else {
                val = (double)(rand() % 20 - 10);
                rowSum += fabs(val);
            }
            A2D[i][j] = val;
            A1D[i * N + j] = val;
        }
        A2D[i][i] += rowSum + 1;
        A1D[i * N + i] += rowSum + 1;
        b[i] = (double)(rand() % 200 - 100);
    }
    printf("Система создана\n\n");

    // ВЫВОД МАТРИЦ (если флаг включён)
    if (SHOW_MATRICES) {
        printf("========== МАТРИЦА A (2D представление) ==========\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%8.2f ", A2D[i][j]);
            }
            printf("| %8.2f\n", b[i]);
        }

        printf("\n========== МАТРИЦА A (1D представление) ==========\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%8.2f ", A1D[i * N + j]);
            }
            printf("| %8.2f\n", b[i]);
        }

        // Проверка, что матрицы одинаковые
        bool identical = true;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (fabs(A2D[i][j] - A1D[i * N + j]) > 1e-10) {
                    identical = false;
                    break;
                }
            }
        }
        printf("\nМатрицы одинаковые: %s\n", identical ? "ДА" : "НЕТ");
        printf("================================================\n\n");
    }

    // Замер для двумерного метода
    printf("Двумерный массив...\n");
    clock_t start = clock();
    x2 = gauss_2d(N, A2D, b);
    clock_t end = clock();
    double time_2d = (double)(end - start) / CLOCKS_PER_SEC;

    // Замер для одномерного метода
    printf("Одномерный массив...\n");
    start = clock();
    x1 = gauss_1d(N, A1D, b);
    end = clock();
    double time_1d = (double)(end - start) / CLOCKS_PER_SEC;



    // Проверка решений (первые 5)
    if (x2 && x1) {
        printf("\nПервые 5 решений (2D):\n");
        for (int i = 0; i < 5 && i < N; i++) {
            printf("x[%d] = %.6f\n", i, x2[i]);
        }

        printf("\nПервые 5 решений (1D):\n");
        for (int i = 0; i < 5 && i < N; i++) {
            printf("x[%d] = %.6f\n", i, x1[i]);
        }

        // Проверка подстановкой для 2D
        printf("\nПроверка (Ax = b) для 2D:\n");
        for (int i = 0; i < 5 && i < N; i++) {
            double sum = 0;
            for (int j = 0; j < N; j++) {
                sum += A2D[i][j] * x2[j];
            }
            printf("Ур%d: %.6f = %.6f\n", i, sum, b[i]);
        }

        // Проверка подстановкой для 1D
        printf("\nПроверка (Ax = b) для 1D:\n");
        for (int i = 0; i < 5 && i < N; i++) {
            double sum = 0;
            for (int j = 0; j < N; j++) {
                sum += A1D[i * N + j] * x1[j];
            }
            printf("Ур%d: %.6f = %.6f\n", i, sum, b[i]);
        }
    }

    // Освобождаем память
    free(A2D);
    free(A1D);
    free(b);
    if (x2) free(x2);
    if (x1) free(x1);


    // Вывод результатов
    printf("\n========================================\n");
    printf("РЕЗУЛЬТАТЫ:\n");
    printf("Двумерный массив: %.4f секунд\n", time_2d);
    printf("Одномерный массив: %.4f секунд\n", time_1d);

    if (time_2d < time_1d) {
        printf("\nДвумерный быстрее на %.2f%%\n", (time_1d - time_2d) / time_1d * 100);
    } else if (time_1d < time_2d) {
        printf("\nОдномерный быстрее на %.2f%%\n", (time_2d - time_1d) / time_2d * 100);
    } else {
        printf("\nВремя одинаковое\n");
    }
    printf("========================================\n");


    return 0;
}