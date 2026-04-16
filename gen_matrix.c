#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char* argv[]) {
    int n = 2000;
    if (argc > 1) n = atoi(argv[1]);
    const char* out = (argc > 2) ? argv[2] : "matrix.txt";

    if (n <= 0) { fprintf(stderr, "❌ Неверный размер\n"); return 1; }

    FILE* f = fopen(out, "w");
    if (!f) { perror("fopen"); return 1; }

    fprintf(f, "%d\n", n);
    srand(42); // Фиксированный seed → одинаковые матрицы при каждом запуске

    for (int i = 0; i < n; i++) {
        double diag_sum = 0.0;
        double *row = (double*)malloc(n * sizeof(double));

        for (int j = 0; j < n; j++) {
            if (i != j) {
                row[j] = (rand() % 201 - 100) / 10.0; // [-10.0; 10.0]
                diag_sum += fabs(row[j]);
            } else {
                row[j] = 0.0; // Заглушка
            }
        }
        // Строгое диагональное преобладание: |a_ii| > Σ|a_ij|
        row[i] = diag_sum + 1.0 + (rand() % 50) / 10.0;
        double b = (rand() % 1001 - 500) / 10.0;

        for (int j = 0; j < n; j++) fprintf(f, "%.6f ", row[j]);
        fprintf(f, "%.6f\n", b);
        free(row);
    }
    fclose(f);
    printf("✅ Сгенерирована матрица %dx%d -> %s\n", n, n, out);
    return 0;
}