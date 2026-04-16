#include "gauss.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 1e-9

void swap_rows_1d(double *matrix, int i, int j, int n_cols, int n) {
    if (i == j) return;
    for (int col = 0; col < n_cols; col++) {
        double temp = matrix[i * n_cols + col];
        matrix[i * n_cols + col] = matrix[j * n_cols + col];
        matrix[j * n_cols + col] = temp;
    }
}

int find_max_row_global_1d(double *local_matrix, int local_rows, int n, int col, int rank, int size, MPI_Comm comm) {
    double local_max_val = 0.0;
    int local_max_row = -1;
    for (int i = 0; i < local_rows; i++) {
        double val = fabs(local_matrix[i * (n + 1) + col]);
        if (val > local_max_val) { local_max_val = val; local_max_row = i; }
    }
    double global_max_val;
    MPI_Allreduce(&local_max_val, &global_max_val, 1, MPI_DOUBLE, MPI_MAX, comm);
    int owner_rank = -1;
    if (local_max_row != -1 && fabs(local_max_val - global_max_val) < EPS) owner_rank = rank;
    int temp_owner = owner_rank;
    MPI_Allreduce(&temp_owner, &owner_rank, 1, MPI_INT, MPI_MAX, comm);
    int global_row = -1;
    if (rank == owner_rank && local_max_row != -1) {
        int rows_per_proc = (n + size - 1) / size;
        global_row = owner_rank * rows_per_proc + local_max_row;
    }
    MPI_Bcast(&global_row, 1, MPI_INT, owner_rank, comm);
    return global_row;
}

void gaussian_elimination_1d(double *local_matrix, int local_rows, int n, int rank, int size, MPI_Comm comm) {
    double *pivot_row = (double *)malloc(n * sizeof(double));
    int rows_per_proc = (n + size - 1) / size;
    int n_cols = n + 1;
    for (int col = 0; col < n; col++) {
        int global_pivot_row = find_max_row_global_1d(local_matrix, local_rows, n, col, rank, size, comm);
        int pivot_owner = global_pivot_row / rows_per_proc;
        int local_pivot_idx = global_pivot_row % rows_per_proc;
        if (rank == pivot_owner && local_pivot_idx < local_rows)
            memcpy(pivot_row, &local_matrix[local_pivot_idx * n_cols], n * sizeof(double));
        MPI_Bcast(pivot_row, n, MPI_DOUBLE, pivot_owner, comm);
        if (fabs(pivot_row[col]) < EPS) { free(pivot_row); MPI_Abort(comm, 1); }
        for (int i = 0; i < local_rows; i++) {
            int global_i = rank * rows_per_proc + i;
            if (global_i == global_pivot_row) continue;
            double current_val = local_matrix[i * n_cols + col];
            if (fabs(current_val) > EPS) {
                double factor = current_val / pivot_row[col];
                for (int j = col; j < n; j++) local_matrix[i * n_cols + j] -= factor * pivot_row[j];
                local_matrix[i * n_cols + col] = 0.0;
            }
        }
    }
    free(pivot_row);
}

void back_substitution_1d(double *local_matrix, int local_rows, int n, int rank, int size, double *solution, MPI_Comm comm) {
    double *x = (double *)malloc(n * sizeof(double));
    int rows_per_proc = (n + size - 1) / size;
    int n_cols = n + 1;
    for (int i = n - 1; i >= 0; i--) {
        int owner = i / rows_per_proc;
        int local_idx = i % rows_per_proc;
        double sum = 0.0;
        if (rank == owner && local_idx < local_rows) {
            for (int j = i + 1; j < n; j++) sum += local_matrix[local_idx * n_cols + j] * x[j];
            x[i] = (local_matrix[local_idx * n_cols + n] - sum) / local_matrix[local_idx * n_cols + i];
        }
        MPI_Bcast(&x[i], 1, MPI_DOUBLE, owner, comm);
    }
    memcpy(solution, x, n * sizeof(double));
    free(x);
}

// ============================================================================
// ГЛАВНАЯ ФУНКЦИЯ (загрузка из файла + распределение)
// ============================================================================
int gauss_1d_parallel(const char* input_file, const char* output_file) {
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 0;
    double *aug_matrix = NULL; // Только для rank 0

    if (rank == 0) {
        double *A = NULL, *b = NULL;
        if (load_system_from_txt(input_file, &n, &A, &b) != 0) { MPI_Abort(MPI_COMM_WORLD, 1); }
        // Формируем расширенную матрицу [A | b]
        aug_matrix = (double*)malloc(n * (n + 1) * sizeof(double));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) aug_matrix[i * (n + 1) + j] = A[i * n + j];
            aug_matrix[i * (n + 1) + n] = b[i];
        }
        free(A); free(b);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int rows_per_proc = (n + size - 1) / size;
    int local_rows = (rank < size - 1) ? rows_per_proc : (n - rank * rows_per_proc);
    if (local_rows < 0) local_rows = 0;

    double *local_matrix = (double *)malloc(local_rows * (n + 1) * sizeof(double));

    // Распределение через MPI_Scatterv
    int *sendcounts = (int*)malloc(size * sizeof(int));
    int *displs = (int*)malloc(size * sizeof(int));
    int offset = 0;
    for (int p = 0; p < size; p++) {
        int r = (p < size - 1) ? rows_per_proc : (n - p * rows_per_proc);
        if (r < 0) r = 0;
        sendcounts[p] = r * (n + 1);
        displs[p] = offset;
        offset += sendcounts[p];
    }
    MPI_Scatterv(aug_matrix, sendcounts, displs, MPI_DOUBLE, local_matrix, local_rows * (n + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) free(aug_matrix);
    free(sendcounts); free(displs);

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    gaussian_elimination_1d(local_matrix, local_rows, n, rank, size, MPI_COMM_WORLD);
    double *solution = (double *)malloc(n * sizeof(double));
    back_substitution_1d(local_matrix, local_rows, n, rank, size, solution, MPI_COMM_WORLD);
    double end = MPI_Wtime();
    double exec_time = end - start;
    double max_time;
    MPI_Reduce(&exec_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        save_solution_to_txt(output_file, n, solution);
        printf("=== Parallel Gaussian (1D) ===\n");
        printf("Size: %d | Procs: %d | Time: %.4f s\n", n, size, max_time);
        printf("Solution saved to: %s\n", output_file);
        // Простая проверка невязки на первых 5 уравнениях
        double max_err = 0.0;
        int check = (local_rows < 5) ? local_rows : 5;
        for (int i = 0; i < check; i++) {
            double comp = 0.0;
            for (int j = 0; j < n; j++) comp += local_matrix[i * (n + 1) + j] * solution[j];
            double err = fabs(comp - local_matrix[i * (n + 1) + n]);
            if (err > max_err) max_err = err;
        }
        printf("Max residual (rank 0): %.2e %s\n", max_err, max_err < 1e-6 ? "✓" : "⚠");
    }

    free(local_matrix); free(solution);
    MPI_Finalize();
    return 0;
}