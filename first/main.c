#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    int rank;
    int count;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count);

    if (rank == 0) {
        char** recv_buffer = malloc(count * sizeof(char*));
        int message_length;
        MPI_Status status;
        for (int i = 0; i < count - 1; i++) {
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_CHAR, &message_length);
            recv_buffer[status.MPI_SOURCE] = malloc(message_length * sizeof(char));
            MPI_Recv(recv_buffer[status.MPI_SOURCE - 1], message_length, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG,
                     MPI_COMM_WORLD, &status);
            printf("%s\n", recv_buffer[status.MPI_SOURCE]);
        }
    } else {
        char str[50];
        sprintf(str, "Hello from %d!", rank);
        MPI_Send(&str, strlen(str) + 1, MPI_CHAR, 0, rank, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
