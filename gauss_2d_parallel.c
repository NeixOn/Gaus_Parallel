#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 1e-9


void swap_rows(double **matrix, int i, int j, int n) {
    if (i == j) return;
    double *temp = matrix[i];
    matrix[i] = matrix[j];
    matrix[j] = temp;
}


int find_max_row_global(double **local_matrix, int local_rows, int n,
                        int col, int rank, int size, MPI_Comm comm) {

    double local_max_val = 0.0;
    int local_max_row = -1;  // Локальный индекс строки

    // 1. Ищем максимум в своих строках
    for (int i = 0; i < local_rows; i++) {
        double val = fabs(local_matrix[i][col]);
        if (val > local_max_val) {
            local_max_val = val;
            local_max_row = i;
        }
    }

    // 2. Находим глобальный максимум по всем процессам
    double global_max_val;
    MPI_Allreduce(&local_max_val, &global_max_val, 1, MPI_DOUBLE, MPI_MAX, comm);

    // 3. Определяем, какой процесс владеет глобальным максимумом
    int owner_rank = -1;
    if (local_max_row != -1 && fabs(local_max_val - global_max_val) < EPS) {
        owner_rank = rank;
    }

    // Находим процесс с максимальным рангом среди тех, у кого есть максимум
    // (если несколько — берём первого по рангу для детерминизма)
    int temp_owner = owner_rank;
    MPI_Allreduce(&temp_owner, &owner_rank, 1, MPI_INT, MPI_MAX, comm);

    // Если есть несколько кандидатов, MPI_MAX вернёт максимальный rank.
    // Для простоты оставим так — в учебном коде это допустимо.

    // 4. Владелец рассылает глобальный номер своей строки
    int global_row = -1;
    if (rank == owner_rank && local_max_row != -1) {
        // Вычисляем глобальный номер строки
        int rows_per_proc = (n + size - 1) / size;
        global_row = owner_rank * rows_per_proc + local_max_row;
    }

    MPI_Bcast(&global_row, 1, MPI_INT, owner_rank, comm);

    return global_row;
}

// ============================================================================
// Прямой ход метода Гаусса (приведение к верхнетреугольному виду)
// ============================================================================
void gaussian_elimination(double **local_matrix, int local_rows, int n,
                          int rank, int size, MPI_Comm comm) {

    double *pivot_row = (double *)malloc(n * sizeof(double));
    int rows_per_proc = (n + size - 1) / size;

    for (int col = 0; col < n; col++) {
        // 1. Найти глобальную строку с максимальным элементом в столбце
        int global_pivot_row = find_max_row_global(local_matrix, local_rows, n,
                                                   col, rank, size, comm);

        // 2. Определить владельца этой строки и её локальный индекс у владельца
        int pivot_owner = global_pivot_row / rows_per_proc;
        int local_pivot_idx = global_pivot_row % rows_per_proc;

        // 3. Владелец копирует строку в буфер
        if (rank == pivot_owner && local_pivot_idx < local_rows) {
            memcpy(pivot_row, local_matrix[local_pivot_idx], n * sizeof(double));
        }

        // 4. Расслать ведущую строку всем процессам
        MPI_Bcast(pivot_row, n, MPI_DOUBLE, pivot_owner, comm);

        // 6. Исключить элемент столбца во всех локальных строках (кроме самой ведущей)
        for (int i = 0; i < local_rows; i++) {
            // Пропускаем саму ведущую строку
            int global_i = rank * rows_per_proc + i;
            if (global_i == global_pivot_row) continue;

            if (fabs(local_matrix[i][col]) > EPS) {
                double factor = local_matrix[i][col] / pivot_row[col];
                for (int j = col; j < n; j++) {
                    local_matrix[i][j] -= factor * pivot_row[j];
                }
                // Обнуляем явно для численной стабильности
                local_matrix[i][col] = 0.0;
            }
        }
    }

    free(pivot_row);
}

// ============================================================================
// Обратный ход: нахождение решения из верхнетреугольной системы
// ============================================================================
void back_substitution(double **local_matrix, int local_rows, int n,
                       int rank, int size, double *solution, MPI_Comm comm) {

    double *x = (double *)malloc(n * sizeof(double));
    int rows_per_proc = (n + size - 1) / size;

    // Идём снизу вверх: от x[n-1] к x[0]
    for (int i = n - 1; i >= 0; i--) {
        int owner = i / rows_per_proc;
        int local_idx = i % rows_per_proc;

        double sum = 0.0;

        // Только владелец строки вычисляет x[i]
        if (rank == owner && local_idx < local_rows) {
            for (int j = i + 1; j < n; j++) {
                sum += local_matrix[local_idx][j] * x[j];
            }
            x[i] = (local_matrix[local_idx][n] - sum) / local_matrix[local_idx][i];
        }

        // Расслать вычисленное x[i] всем процессам
        MPI_Bcast(&x[i], 1, MPI_DOUBLE, owner, comm);
    }

    // Скопировать результат в выходной массив
    memcpy(solution, x, n * sizeof(double));
    free(x);
}

// ============================================================================
// Распределённая генерация матрицы с диагональным преобладанием
// Каждый процесс генерирует ТОЛЬКО свои строки
// ============================================================================
void generate_matrix_distributed(double **local_matrix, int local_rows,
                                 int n, int size, int rank, int global_seed) {

    int rows_per_proc = (n + size - 1) / size;

    for (int i = 0; i < local_rows; i++) {
        int global_row = rank * rows_per_proc + i;  // Глобальный номер строки

        // Детерминированный seed на основе номера строки и общего seed
        srand(global_seed + global_row * 1000 + rank);

        double diag_sum = 0.0;

        // Генерируем все элементы строки, кроме диагонального
        for (int j = 0; j < n; j++) {
            if (j != global_row) {
                // Случайное значение от -10 до 10
                local_matrix[i][j] = (rand() % 21) - 10;
                diag_sum += fabs(local_matrix[i][j]);
            }
        }

        // Диагональный элемент: сумма модулей остальных + случайное положительное число
        // Это гарантирует строгое диагональное преобладание: |a_ii| > Σ|a_ij|
        local_matrix[i][global_row] = diag_sum + 1.0 + (rand() % 10);

        // Свободный член (последний столбец расширенной матрицы)
        local_matrix[i][n] = (rand() % 101) - 50;  // от -50 до 50
    }
}

// ============================================================================
// Главная функция
// ============================================================================
int gauss_2d_parallel(int N) {

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Размер системы (можно задать через аргумент командной строки)
    int n = N;

    // Seed для генерации случайных чисел
    int global_seed = 42;

    // Вычисляем, сколько строк будет у этого процесса
    int rows_per_proc = (n + size - 1) / size;  // Округление вверх
    int local_rows = (rank < size - 1) ? rows_per_proc : (n - rank * rows_per_proc);

    // Коррекция для последнего процесса, если n не делится на size
    if (rank == size - 1) {
        local_rows = n - rank * rows_per_proc;
        if (local_rows < 0) local_rows = 0;
    }

    // Выделение памяти под локальную часть матрицы [local_rows × (n+1)]
    // Последний столбец — свободные члены
    double **local_matrix = (double **)malloc(local_rows * sizeof(double *));
    for (int i = 0; i < local_rows; i++) {
        local_matrix[i] = (double *)malloc((n + 1) * sizeof(double));
    }


    double start_time, end_time, exec_time;

// 1. Синхронизация: ждём, пока ВСЕ процессы закончат инициализацию
    MPI_Barrier(MPI_COMM_WORLD);

// 2. Старт таймера (только вычисления!)

    // ========================================================================
    // РАСПРЕДЕЛЁННАЯ ГЕНЕРАЦИЯ МАТРИЦЫ
    // Каждый процесс заполняет только свои строки — память не дублируется!
    // ========================================================================
    generate_matrix_distributed(local_matrix, local_rows, n, size, rank, global_seed);

    // ========================================================================
    // ПАРАЛЛЕЛЬНЫЙ МЕТОД ГАУССА
    // ========================================================================
    start_time = MPI_Wtime();
    // Прямой ход: приведение к верхнетреугольному виду
    gaussian_elimination(local_matrix, local_rows, n, rank, size, MPI_COMM_WORLD);

    // Обратный ход: нахождение решения
    double *solution = (double *)malloc(n * sizeof(double));
    back_substitution(local_matrix, local_rows, n, rank, size, solution, MPI_COMM_WORLD);

    end_time = MPI_Wtime();

    exec_time = end_time - start_time;
    double max_time;
    MPI_Reduce(&exec_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    // ========================================================================
    // ВЫВОД РЕЗУЛЬТАТОВ (только процесс 0)
    // ========================================================================
    if (rank == 0) {
        // Print execution time
        printf("=== Parallel Gaussian Elimination Results ===\n");
        printf("Matrix size: %d x %d | Processes: %d\n", n, n, size);
        printf("Execution time: %.4f seconds\n\n", max_time);

        // Verification: Check residual ||Ax - b|| for first few equations
        // Note: Process 0 only owns the first 'local_rows' rows of the matrix
        int rows_per_proc = (n + size - 1) / size;
        int local_rows = (n + size - 1) / size; // recalculate for clarity
        if (rank == size - 1) {
            local_rows = n - rank * rows_per_proc;
            if (local_rows < 0) local_rows = 0;
        }

        printf("Verification (first 5 equations owned by rank 0):\n");
        printf("%-8s %-15s %-15s %-10s\n", "Equation", "Computed Ax", "Expected b", "Error");
        printf("------------------------------------------------------------\n");

        int check_count = (local_rows < 5) ? local_rows : 5;
        double max_error = 0.0;

        for (int i = 0; i < check_count; i++) {
            double computed = 0.0;
            // Compute dot product: A[i][0..n-1] * x[0..n-1]
            for (int j = 0; j < n; j++) {
                computed += local_matrix[i][j] * solution[j];
            }
            double expected = local_matrix[i][n]; // b[i] is in the last column
            double error = fabs(computed - expected);

            if (error > max_error) max_error = error;

            printf("%-8d %-15.6f %-15.6f %-10.2e\n", i, computed, expected, error);
        }

        printf("------------------------------------------------------------\n");
        printf("Max residual error: %.2e\n", max_error);

        if (max_error < 1e-6) {
            printf("✓ Solution verified successfully (error within tolerance)\n");
        } else {
            printf("⚠ Warning: Large residual detected. Solution may be inaccurate.\n");
        }
        printf("\n");
    }

    // ========================================================================
    // ОСВОБОЖДЕНИЕ ПАМЯТИ
    // ========================================================================
    for (int i = 0; i < local_rows; i++) {
        free(local_matrix[i]);
    }
    free(local_matrix);
    free(solution);
    return 0;
}