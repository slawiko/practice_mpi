#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>

double f(double x) {
    return x + 1;
}

double calculate_block(double* intervals, int count, double interval_size) {
    double result = 0.0, a, b;

    for (int i = 0; i < count; ++i) {
        a = intervals[i] - interval_size;
        b = intervals[i];
        result += (f(a) + f(b)) / 2. * (b - a);
    }

    return result;
}

int main(int argc, char* argv[]) {
    int a = 0, b = 10;
    int tag_ready = 1;
    int tag_new_block = 2;
    int tag_result = 3;
    int tag_finish = 4;

    MPI_Init(&argc, &argv);
    int rank;
    int p;

    int n, r, block_count;
    n = atoi(argv[1]);
    r = atoi(argv[2]);

    block_count = n / r;
    double interval_size = (b - a) / n;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (rank == 0) {
        double* intervals = malloc(sizeof(double) * n);

        intervals[0] = interval_size;
        for (int i = 1; i < n; ++i) {
            intervals[i] = intervals[i - 1] + interval_size;
        }

        MPI_Status status;


        for (int i = 0; i < block_count; ++i) {
            MPI_Recv(&tag_ready, 1, MPI_INT, MPI_ANY_SOURCE, tag_ready, MPI_COMM_WORLD, &status);
            MPI_Send(intervals + i * r, r, MPI_DOUBLE, status.MPI_SOURCE, tag_new_block, MPI_COMM_WORLD);
        }

        double temp = 0.0;
        for (int i = 0; i < p - 1; ++i) {
            MPI_Recv(&tag_ready, 1, MPI_INT, MPI_ANY_SOURCE, tag_ready, MPI_COMM_WORLD, &status);
            MPI_Send(&temp, 1, MPI_DOUBLE, status.MPI_SOURCE, tag_finish, MPI_COMM_WORLD);
        }

        double result = 0.0, result_buffer;
        for (int i = 0; i < p - 1; ++i) {
            MPI_Recv(&result_buffer, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag_result, MPI_COMM_WORLD, &status);
            result += result_buffer;
        }

        free(intervals);
        printf("%f", result);
    } else {
        MPI_Status status;

        double result = 0.0;
        double* block = malloc(sizeof(double) * r);
        while (true) {
            MPI_Send(&tag_ready, 1, MPI_INT, 0, tag_ready, MPI_COMM_WORLD);
            MPI_Recv(block, r, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == tag_new_block) {
                result += calculate_block(block, r, interval_size);
                continue;
            }
            if (status.MPI_TAG == tag_finish) {
                MPI_Send(&result, 1, MPI_DOUBLE, 0, tag_result, MPI_COMM_WORLD);
                break;
            }
        }

        free(block);
    }

    MPI_Finalize();

    return 0;
}