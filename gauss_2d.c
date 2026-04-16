#include "gauss.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double* gauss_2d(int n, double **A, double *b) {
    double *x = malloc(n * sizeof(double));

    // ✅ Выделяем память как массив указателей (вместо VLA)
    double **a = (double **)calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++) {
        a[i] = (double *)calloc(n, sizeof(double));
    }
    double *bb = malloc(n * sizeof(double));

    // Копируем данные
    for (int i = 0; i < n; i++) {
        bb[i] = b[i];
        for (int j = 0; j < n; j++) {
            a[i][j] = A[i][j];
        }
    }

    // Прямой ход
    for (int k = 0; k < n - 1; k++) {
        // Поиск главного элемента
        int maxRow = k;
        double maxVal = fabs(a[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(a[i][k]) > maxVal) {
                maxVal = fabs(a[i][k]);
                maxRow = i;
            }
        }

        // Перестановка строк
        if (maxRow != k) {
            for (int j = k; j < n; j++) {
                double tmp = a[k][j];
                a[k][j] = a[maxRow][j];
                a[maxRow][j] = tmp;
            }
            double tmp = bb[k];
            bb[k] = bb[maxRow];
            bb[maxRow] = tmp;
        }

        // Проверка на вырожденность
        if (fabs(a[k][k]) < 1e-12) {
            printf("Ошибка: матрица вырождена на шаге %d\n", k);
            for (int i = 0; i < n; i++) free(a[i]);
            free(a);
            free(bb);
            free(x);
            return NULL;
        }

        // Обнуление
        for (int i = k + 1; i < n; i++) {
            double factor = a[i][k] / a[k][k];
            for (int j = k; j < n; j++) {
                a[i][j] -= factor * a[k][j];
            }
            bb[i] -= factor * bb[k];
        }
    }

    // Проверка последней диагонали
    if (fabs(a[n-1][n-1]) < 1e-12) {
        printf("Ошибка: матрица вырождена!\n");
        for (int i = 0; i < n; i++) free(a[i]);
        free(a);
        free(bb);
        free(x);
        return NULL;
    }

    // Обратный ход
    for (int i = n - 1; i >= 0; i--) {
        double sum = bb[i];
        for (int j = i + 1; j < n; j++) {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }

    // ✅ Освобождаем память (двумерный массив)
    for (int i = 0; i < n; i++) free(a[i]);
    free(a);
    free(bb);
    return x;
}