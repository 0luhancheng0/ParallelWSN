#include <stdio.h>
#include <memory.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include <unistd.h>
#include "wsn.h"
#include "cipher.h"
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
    int global_size, global_rank, size, rank, nneighbor;
    uint8_t random_num;
    uint8_t* message;
    MPI_Comm node_comm;
    int r0, r1;
    const int baserank_list[1] = {BASERANK};
    int dim[2] = {X_SIZE, Y_SIZE}, period[2] = {0}, coord[2], reorder = 0;
    int* neighbor_rank;
//	MPI_Datatype message_type;
	
    MPI_Group global_group, node_group;
    MPI_Request *reqs;
    MPI_Request trash_req;

    MPI_Init(&argc, &argv);
//	MPI_Type_contiguous(MESSAGE_LEN, MPI_UINT8_T, &message_type);
//	MPI_Type_commit(&message_type);
    uint8_t* neighbor_result;
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


        // create group for inter node communication
        MPI_Comm_create_group(MPI_COMM_WORLD, node_group, 0, &node_comm);
        MPI_Group_rank(node_group, &rank);
        MPI_Group_size(node_group, &size);

        // create cart topology
        MPI_Cart_create(node_comm, 2, dim, period, reorder, &node_comm);
        MPI_Cart_coords(node_comm, rank, 2, coord);
        get_neighbor_count(coord, &nneighbor);

        // setup different random number seed on different rank
        srand(time(NULL)+rank);
        random_num = (uint8_t)(rand() & ((1 << N_BIT_RAND) - 1));
        printf("Random number is %d on rank %d coordinates are %d %d. Has %d neighbor\n", random_num, rank, coord[0], coord[1], nneighbor);

        // add padding of zeros following random num
        message = (uint8_t*)calloc(MESSAGE_LEN, 1);
        memcpy(message, &random_num, 1);

        encrypt(message, message, MESSAGE_LEN, key, KEY_SIZE);
        

        MPI_Alloc_mem((MPI_Aint)(nneighbor * MESSAGE_LEN), MPI_INFO_NULL, &neighbor_result);
        MPI_Alloc_mem((MPI_Aint)(sizeof(MPI_Request) * 4), MPI_INFO_NULL, &reqs);
        MPI_Alloc_mem(nneighbor, MPI_INFO_NULL, &neighbor_rank);
		
		// similar implementation with MPI_Ineighbor_alltoall
        for (int i=0, dim=0; dim<2; ++dim) {
            MPI_Cart_shift(node_comm, dim, 1, &r0, &r1);
            if (r0 >= 0) {
                MPI_Isend(message, MESSAGE_LEN, MPI_UINT8_T, r0, 0, node_comm, &trash_req);
                MPI_Irecv(&neighbor_result[i*MESSAGE_LEN], MESSAGE_LEN, MPI_UINT8_T, r0, MPI_ANY_TAG, node_comm, &reqs[i]);
                neighbor_rank[i] = r0;
                i++;
            }
            if (r1 >= 0) {
                MPI_Isend(message, MESSAGE_LEN, MPI_UINT8_T, r1, 0, node_comm, &trash_req);
                MPI_Irecv(&neighbor_result[i*MESSAGE_LEN], MESSAGE_LEN, MPI_UINT8_T, r1, MPI_ANY_TAG, node_comm, &reqs[i]);
                neighbor_rank[i] = r1;
                i++;
            }
        }
        MPI_Waitall(nneighbor, reqs, MPI_STATUSES_IGNORE);
        decrypt(neighbor_result, neighbor_result, nneighbor * MESSAGE_LEN, key, KEY_SIZE);
        
        for (int i=0;i <nneighbor; i++) {
            printf("rank %d received %d from %d\n", rank, neighbor_result[i*MESSAGE_LEN], neighbor_rank[i]);
        }
		MPI_Comm_free(&node_comm);
        free(message);
    } else {
		printf("Rank of base station is %d\n", global_rank);
	}
    
    MPI_Finalize();
    return(0);
}
