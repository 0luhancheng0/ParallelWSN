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

// The time taken to encrypt/decrypt a message before sending or after receiving. The use of OpenMP here should demonstrate improvements in encryption and/or decryption processing time.
// - Number of messages passed throughout the network
// - Number of events occurred throughout the network
// - Details of nodes involved in each of the events (reference node and its adjacent nodes)

struct event {
	int occur_on_ranks[4];
    double encryption_time;
	double decryption_time;
    long int timestamp;
    int n_times;
    int iteration;
    int reference_rank;
	uint8_t num;
    
};

struct logger {
	int encrypt_t;
	int decryp_t;
	int n_message;
	int n_event;
	int* reference_node;
	int** adjacent_node;
};

void print_event(struct event e) {
    printf("event num %d detected on rank (local) %d\nencryption time taken: %lf seconds\ndecryption time taken: %lf seconds\ntime: %siteration: %d\n",e.num, e.reference_rank, e.encryption_time, e.decryption_time, asctime(localtime(&e.timestamp)), e.iteration);
    printf("number generated %d times on rank: ", e.n_times);
    for (int i=0;i<e.n_times;i++) {
        printf("%d ", e.occur_on_ranks[i]);
    }
    printf("\n\n");
    
}


static void get_neighbor_count(int coord[2], int* neighbor_count) {
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
    struct event e;
    // create mpi struct type for event
    // http://www.catb.org/esr/structure-packing/
    int event_property_count = 7;
    int blocklengths[] = {4, 1, 1, 1, 1, 1, 1, 1};
    
    MPI_Aint displacement[] = {0, 16, 24, 32, 40, 44, 48, 52};
    MPI_Datatype types[] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_INT, MPI_INT, MPI_INT, MPI_UINT8_T};
    MPI_Datatype my_mpi_event_type;

    
    struct event *event_recv_buff;

    // MPI_Type_struct(8, );
    
    double tick, encryption_time, decryption_time, t0;
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
    MPI_Type_create_struct(event_property_count, blocklengths, displacement, types, &my_mpi_event_type);
    MPI_Type_commit(&my_mpi_event_type);
    
    tick = MPI_Wtick();
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

        message = calloc(MESSAGE_LEN, 1);
        MPI_Alloc_mem((MPI_Aint)(nneighbor * MESSAGE_LEN), MPI_INFO_NULL, &neighbor_result);
        MPI_Alloc_mem((MPI_Aint)(sizeof(MPI_Request) * 4), MPI_INFO_NULL, &reqs);
        MPI_Alloc_mem(nneighbor, MPI_INFO_NULL, &neighbor_rank);

        // same as power(2, nbit), but more efficient
        int upperbound = 1 << N_BIT_RAND;
        int *num_recv = calloc(upperbound, sizeof(int));
        // e.occur_on_ranks = malloc(nneighbor);

        for (int current_i = 0; current_i < N_ITERATION; current_i++)
        {

            random_num = (uint8_t)(rand() & ((1 << N_BIT_RAND) - 1));
            printf("Random number is %d on rank %d coordinates are %d %d. Has %d neighbor\n", random_num, rank, coord[0], coord[1], nneighbor);

            // add padding of zeros following random num
            
            memcpy(message, &random_num, 1);

            t0 = MPI_Wtime(); 
            encrypt(message, message, MESSAGE_LEN, key, KEY_SIZE);
            encryption_time = MPI_Wtime() - t0;

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


            t0 = MPI_Wtime();
            decrypt(neighbor_result, neighbor_result, nneighbor * MESSAGE_LEN, key, KEY_SIZE);
            decryption_time = MPI_Wtime() - t0;


            // int event_occured = 0;
            
            for (int i=0;i <nneighbor; i++) {
                // printf("rank %d received %d from %d\n", rank, neighbor_result[i*MESSAGE_LEN], neighbor_rank[i]);
                num_recv[(int)neighbor_result[i*MESSAGE_LEN]] ++; 
            }
            // MPI_Barrier(node_comm);
            for (int i=0;i < upperbound; i++) {
                // printf("num_recv[i]=%d on rank %d\n", i, num_recv[i], rank);
                if (num_recv[i] >= RAND_TH)
                {
                    e.num = (uint8_t)i;
                    e.n_times = num_recv[i];
                    e.iteration = current_i;
                    e.reference_rank = rank;
                    e.encryption_time = encryption_time;
                    e.decryption_time = decryption_time;
                    e.timestamp = time(NULL);

                    
                    // int k=0;
                    // printf("nneighbor = %d\n", nneighbor);
                    for (int j = 0, k=0; j < nneighbor; j++)
                    {
                        // printf("k = %d, j = %d, neighbor_result[%d * MESSAGE_LEN] = %d\n", k, j, i, neighbor_result[i*MESSAGE_LEN]);
                        if ((int)neighbor_result[j * MESSAGE_LEN] == i)
                        {
                            e.occur_on_ranks[k++] = neighbor_rank[j];
                            // printf("k=%d\n", k);
                        }
                    }
                    // print_event(e);
                    MPI_Send(&e, 1, my_mpi_event_type, BASERANK, 0, MPI_COMM_WORLD);
                }

                // clean up for next iteration
                num_recv[i] = 0;
            }
            usleep(1000 * INTERVAL);
        }
		MPI_Comm_free(&node_comm);
        free(message);
    } else {
        printf("bask rank is %d\n", global_rank);
        event_recv_buff = malloc(sizeof(struct event) * X_SIZE * Y_SIZE);


        MPI_Recv(event_recv_buff, 1, my_mpi_event_type, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        print_event(*event_recv_buff);
        // printf("printed by base rank\n");
    }
    MPI_Finalize();
    return(0);
    
}


