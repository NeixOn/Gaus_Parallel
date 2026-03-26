#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include "gauss.h"

// Создаёт случайную матрицу с диагональным преобладанием
// Возвращает указатель на матрицу (2D) и правую часть
void create_system(int n, double (*A2D)[n], double *A1D, double *b) {
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
    // Выделяем память
    double (*A2D)[n] = malloc(sizeof(double[n][n]));
    double *A1D = malloc(n * n * sizeof(double));
    double *b = malloc(n * sizeof(double));
    double *x2 = NULL;
    double *x1 = NULL;

    if (!A2D || !A1D || !b) {
        //printf("Ошибка выделения памяти\n");
        free(A2D);
        free(A1D);
        free(b);
        return;
    }

    //printf("Создание системы %d x %d...\n", n, n);
    srand(time(NULL));
    create_system(n, A2D, A1D, b);
    //printf("Система создана\n\n");

    // Вывод матриц
    if (show_matrices) {
        //printf("========== МАТРИЦА A (2D представление) ==========\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%8.2f ", A2D[i][j]);
            }
            printf("| %8.2f\n", b[i]);
        }

        //printf("\n========== МАТРИЦА A (1D представление) ==========\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                //printf("%8.2f ", A1D[i * n + j]);
            }
            //printf("| %8.2f\n", b[i]);
        }
        //printf("\n");
    }

    // Замер для 2D метода
    //printf("Двумерный массив...\n");
    clock_t start = clock();
    x2 = gauss_2d(n, A2D, b);
    clock_t end = clock();
    double time_2d = (double)(end - start) / CLOCKS_PER_SEC;

    // Замер для 1D метода
    //printf("Одномерный массив...\n");
    start = clock();
    x1 = gauss_1d(n, A1D, b);
    end = clock();
    double time_1d = (double)(end - start) / CLOCKS_PER_SEC;



//    if (time_2d < time_1d) {
//        printf("\nДвумерный быстрее на %.2f%%\n", (time_1d - time_2d) / time_1d * 100);
//    } else if (time_1d < time_2d) {
//        printf("\nОдномерный быстрее на %.2f%%\n", (time_2d - time_1d) / time_2d * 100);
//    } else {
//        printf("\nВремя одинаковое\n");
//    }
//    printf("========================================\n");

//    // Проверка решений (первые 5)
//    if (x2 && x1) {
//        printf("\nПервые 5 решений (2D):\n");
//        for (int i = 0; i < 5 && i < n; i++) {
//            printf("x[%d] = %.6f\n", i, x2[i]);
//        }
//
//        printf("\nПервые 5 решений (1D):\n");
//        for (int i = 0; i < 5 && i < n; i++) {
//            printf("x[%d] = %.6f\n", i, x1[i]);
//        }
//
//        // Проверка подстановкой
//        printf("\nПроверка (Ax = b) для 2D:\n");
//        for (int i = 0; i < 5 && i < n; i++) {
//            double sum = 0;
//            for (int j = 0; j < n; j++) {
//                sum += A2D[i][j] * x2[j];
//            }
//            printf("Ур%d: %.6f = %.6f\n", i, sum, b[i]);
//        }
//
//        printf("\nПроверка (Ax = b) для 1D:\n");
//        for (int i = 0; i < 5 && i < n; i++) {
//            double sum = 0;
//            for (int j = 0; j < n; j++) {
//                sum += A1D[i * n + j] * x1[j];
//            }
//            printf("Ур%d: %.6f = %.6f\n", i, sum, b[i]);
//        }
//    }

    // Освобождаем память
    free(A2D);
    free(A1D);
    free(b);
    if (x2) free(x2);
    if (x1) free(x1);


    // Вывод результатов
    printf("\n========================================\n");
    printf("Results:\n");
    printf("2d: %.4f second\n", time_2d);
    printf("1d: %.4f second\n", time_1d);
}