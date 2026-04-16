#include "gauss.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]) {
    const char* input  = (argc > 1) ? argv[1] : "matrix.txt";
    const char* output = (argc > 2) ? argv[2] : "solution_seq.txt";

    int n; double *A = NULL, *b = NULL;
    printf("📥 Загрузка: %s\n", input);
    if (load_system_from_txt(input, &n, &A, &b) != 0) return 1;

    // 1D версия
    double *x1 = gauss_1d(n, A, b);
    if (x1) {
        save_solution_to_txt(output, n, x1);
        printf("✅ 1D решение сохранено в %s\n", output);
        free(x1);
    }

    // 2D версия
    double **A2d = flat_to_2d(n, A);
    double *x2 = gauss_2d(n, A2d, b);
    if (x2) {
        char out2[256]; snprintf(out2, sizeof(out2), "%s_2d.txt", output);
        save_solution_to_txt(out2, n, x2);
        printf("✅ 2D решение сохранено в %s\n", out2);
        free(x2);
    }
    free(A2d); free(A); free(b);
    return 0;
}