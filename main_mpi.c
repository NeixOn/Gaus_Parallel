#include <stdio.h>
#include <stdlib.h>
#include "gauss.h"

int main(int argc, char* argv[]) {
    const char* in  = (argc > 1) ? argv[1] : "matrix.txt";
    const char* out1 = (argc > 2) ? argv[2] : "solution_1d_mpi.txt";
    const char* out2 = (argc > 3) ? argv[3] : "solution_2d_mpi.txt";

    printf("🚀 Запуск параллельных версий...\n");
    gauss_1d_parallel(in, out1);
    // Запуск 2D версии в том же пуле процессов требует отдельного вызова или реинициализации MPI.
    // Для простоты запускайте их отдельно или раскомментируйте нужную версию.
    // gauss_2d_parallel(in, out2); 
    return 0;
}