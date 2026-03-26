// В main.c:
#include <mpi.h>
#include "gauss.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);  // ✅ Инициализация в главном
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int N = 2000;

    // Запустите только нужный алгоритм (или оба, если нужно)
    gauss_2d_parallel(N);  // ✅ Без MPI_Init внутри
    gauss_1d_parallel(N); // ✅ Без MPI_Init внутри

    // Если run_test тоже требует MPI — запустите её здесь
    // run_test(N, true);  // Если нужна параллельная версия
    if (rank == 0) {
        run_test(N, false);
    }
    MPI_Finalize();
    return 0;
}