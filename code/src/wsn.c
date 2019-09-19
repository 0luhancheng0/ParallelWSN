#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>
#include "config.h"
typedef struct logging
{
    int n_event;
    int* enent_num;
    int* reference_rank;
} logging;

struct event_message 
{
    int event_num;
    int my_rank;
    int* reference_rank;
};
void get_neighbor_count(int coord[2], int* neighbor_count) {
	int *neibhbor_count;
	// MPI_Cart_shift()	
    *neighbor_count = 4;
    if (coord[0] == 0) {
        *neighbor_count -= 1;
    }
    if (coord[0] == X_SIZE-1) {
        *neighbor_count -= 1;
    }

    if (coord[1] == 0) {
        *neighbor_count -= 1;
    }
    if (coord[1] == Y_SIZE - 1)
    {
        *neighbor_count -= 1;
    }
}


int main(int argc, char *argv[])
{
    int global_size, global_rank, size, rank, nneighbor, colored_rank, colored_size;
    uint8_t random_num;
    uint8_t* neighbor_result;
    MPI_Comm node_comm;
    int r0, r1;
    const int baserank_list[1] = {BASERANK};
    int dim[2] = {X_SIZE, Y_SIZE}, period[2] = {0}, coord[2], reorder = 0;
    int* neighbor_rank;
    MPI_Group global_group, node_group;
    MPI_Request *reqs;
    // MPI_Request reqs[4];
    MPI_Request trash_req;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    
    if (global_size != X_SIZE * Y_SIZE + 1)
    {
        printf("Please run with %d processes.\n", X_SIZE * Y_SIZE + 1);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Comm_group(MPI_COMM_WORLD, &global_group);
    MPI_Group_excl(global_group, 1, baserank_list, &node_group);

    if (global_rank != BASERANK) {
        MPI_Comm_create_group(MPI_COMM_WORLD, node_group, 0, &node_comm);
        MPI_Group_rank(node_group, &rank);
        MPI_Group_size(node_group, &size);
        MPI_Cart_create(node_comm, 2, dim, period, reorder, &node_comm);
        MPI_Cart_coords(node_comm, rank, 2, coord);
        srand(time(NULL)+rank);

        random_num = (uint8_t)(rand() & ((1 << N_BIT_RAND) - 1));

        get_neighbor_count(coord, &nneighbor);

        printf("Random number is %d on rank %d coordinates are %d %d. Has %d neighbor\n", random_num, rank, coord[0], coord[1], nneighbor);
        MPI_Alloc_mem((MPI_Aint)(nneighbor * sizeof(MPI_UINT8_T)), MPI_INFO_NULL, &neighbor_result);
        MPI_Alloc_mem((MPI_Aint)(sizeof(MPI_Request) * 4), MPI_INFO_NULL, &reqs);
        MPI_Alloc_mem(nneighbor, MPI_INFO_NULL, &neighbor_rank);
        for (int i=0, dim=0; dim<2; ++dim) {
            MPI_Cart_shift(node_comm, dim, 1, &r0, &r1);
            if (r0 >= 0) {
                MPI_Isend(&random_num, 1, MPI_UINT8_T, r0, 0, node_comm, &trash_req);
                MPI_Irecv(&neighbor_result[i], 1, MPI_UINT8_T, r0, MPI_ANY_TAG, node_comm, &reqs[i]);
                neighbor_rank[i] = r0;
                i++;
            }
            if (r1 >= 0) {
                MPI_Isend(&random_num, 1, MPI_UINT8_T, r1, 0, node_comm, &trash_req);
                MPI_Irecv(&neighbor_result[i], 1, MPI_UINT8_T, r1, MPI_ANY_TAG, node_comm, &reqs[i]);
                neighbor_rank[i] = r1;
                i++;
            }
        }
        MPI_Waitall(nneighbor, reqs, MPI_STATUSES_IGNORE);
        for (int i=0;i <nneighbor; i++) {
            printf("rank %d received %u from %d\n", rank, neighbor_result[i], neighbor_rank[i]);
        }
		// MPI_Free_mem(&neighbor_result);
		// MPI_Free_mem(&neighbor_rank);
		// MPI_Free_mem(&reqs);
		MPI_Comm_free(&node_comm);
        // MPI_Finalize();
        // return(0);
    }
    
    MPI_Finalize();
    return(0);
}
