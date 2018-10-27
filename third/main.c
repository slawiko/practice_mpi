#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int first_row_tag = 0;

    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int r = 5, N1 = 25, N2 = 25;
    int block_cnt = N2 / r;
    int block_i = N1 / p;

    int cur_block = 0;

    int **block = malloc(sizeof(int *) * block_i);
    for (int i = 0; i < block_i; ++i) {
        block[i] = malloc(sizeof(int) * r);
    }

    MPI_Status status;

    int *first_row = malloc(sizeof(int) * r);
    int *first_row_2 = malloc(sizeof(int) * r);
    MPI_Request first_row_request;

    int *last_row = malloc(sizeof(int) * r);
    MPI_Request last_row_request;
    if (rank != 0) {
        MPI_Irecv(first_row, r, MPI_INT, rank - 1, first_row_tag, MPI_COMM_WORLD, &first_row_request);
    }

    while (cur_block < block_cnt) {
        if (rank == 0) {
            // init first_row with 0
            for (int j = 0; j < r; ++j) {
                block[0][j] = 0;
            }
        } else {
            // recv first_row
            MPI_Wait(&first_row_request, &status);
            for (int j = 0; j < r; ++j) {
                block[0][j] = first_row[j];
            }
        }

        if (cur_block != block_cnt - 1) {
            MPI_Irecv(first_row, r, MPI_INT, rank - 1, first_row_tag, MPI_COMM_WORLD, &first_row_request);
        }
        // calculate block
        for (int i = 1; i < block_i; ++i) {
            for (int j = 0; j < r; ++j) {
                block[i][j] = block[i - 1][j] + 1;
            }
        }

        printf("I'm %d, block %d\n", rank, cur_block);
        for (int i = 0; i < block_i; ++i) {
            for (int j = 0; j < r; ++j) {
                printf("%d ", block[i][j]);
            }
            printf("\n");
        }

        if (rank + 1 != p) {
            if (cur_block != 0) {
                MPI_Wait(&last_row_request, &status);
            }

            // copy last row
            for (int j = 0; j < r; ++j) {
                last_row[j] = block[r - 1][j];
            }

            // send last_row
            MPI_Isend(last_row, r, MPI_INT, rank + 1, first_row_tag, MPI_COMM_WORLD, &last_row_request);
        }

        cur_block++;
    }


    MPI_Finalize();

    return 0;
}