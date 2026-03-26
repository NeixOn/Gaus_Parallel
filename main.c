// В main.c:
#include <mpi.h>
#include <stdlib.h>
#include "gauss.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int N = 4000;
    int mode = 0;  // 0=2d, 1=1d, 2=seq

    if (argc > 1) N = atoi(argv[1]);
    if (argc > 2) mode = atoi(argv[2]);

    if (mode == 0) {
        gauss_2d_parallel(N);
    } else if (mode == 1) {
        gauss_1d_parallel(N);
    } else if (mode == 2 && rank == 0) {
        run_test(N, false);
    }

    MPI_Finalize();
    return 0;
}