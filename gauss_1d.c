#include "gauss.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double* gauss_1d(int n, double *A, double *b){
    double *x = malloc(n * sizeof(double));
    // Выделяем память под копии
    double *a = malloc(n * n * sizeof(double));
    double *bb = malloc(n * sizeof(double));

    // Копируем данные
    for (int i = 0; i < n; i++) {
        bb[i] = b[i];
        for (int j = 0; j < n; j++) {
            a[i * n + j] = A[i * n + j];
        }
    }

    // Прямой ход
    for (int k = 0; k < n - 1; k++) {
        // Поиск главного элемента
        int maxRow = k;
        double maxVal = fabs(a[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(a[i * n + k]) > maxVal) {
                maxVal = fabs(a[i * n + k]);
                maxRow = i;
            }
        }

        // Перестановка строк
        if (maxRow != k) {
            for (int j = k; j < n; j++) {
                double tmp = a[k * n + j];
                a[k * n + j] = a[maxRow * n + j];
                a[maxRow * n + j] = tmp;
            }
            double tmp = bb[k];
            bb[k] = bb[maxRow];
            bb[maxRow] = tmp;
        }

        // Проверка на вырожденность
        if (fabs(a[k * n + k]) < 1e-12) {
            printf("Ошибка: матрица вырождена на шаге %d\n", k);
            free(a);
            free(bb);
            return 0;
        }

        // Исключение
        for (int i = k + 1; i < n; i++) {
            double factor = a[i * n + k] / a[k * n + k];
            for (int j = k; j < n; j++) {
                a[i * n + j] -= factor * a[k * n + j];
            }
            bb[i] -= factor * bb[k];
        }
    }

    // Проверка последней диагонали
    if (fabs(a[(n-1) * n + (n-1)]) < 1e-12) {
        printf("Ошибка: матрица вырождена!\n");
        free(a);
        free(bb);
        return 0;
    }

    // Обратный ход
    for (int i = n - 1; i >= 0; i--) {
        double sum = bb[i];
        for (int j = i + 1; j < n; j++) {
            sum -= a[i * n + j] * x[j];
        }
        x[i] = sum / a[i * n + i];
    }

    free(a);
    free(bb);
    return x;
}