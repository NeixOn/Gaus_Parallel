#ifndef GAUSS_H
#define GAUSS_H
#include <stdbool.h>

// Последовательные версии
double* gauss_2d(int n, double **A, double *b);
double* gauss_1d(int n, double *A, double *b);

// Параллельные версии (принимают имена файлов)
int gauss_2d_parallel(const char* input_file, const char* output_file);
int gauss_1d_parallel(const char* input_file, const char* output_file);

// Файловый ввод/вывод
int load_system_from_txt(const char* filename, int* n, double** A, double** b);
int save_solution_to_txt(const char* filename, int n, const double* x);

// Утилиты
double** flat_to_2d(int n, double* A_flat);

// Тесты
void run_test(int n, bool show_matrices);

#endif