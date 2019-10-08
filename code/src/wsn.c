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
#include <assert.h>

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
	int num;
    int reference_rank;
};

struct logger {
	int encrypt_time;
	int decryp_time;
	int n_message;
	int n_event;
	int* reference_node;
	int** adjacent_node;
};

static void print_event(struct event e) {
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
// 
// struct event {
// 	int occur_on_ranks[4];
//     double encryption_time;
// 	double decryption_time;
//     long int timestamp;
//     int n_times;
//     int iteration;
//     int reference_rank;
// 	uint8_t num;
//     
// };

int main(int argc, char *argv[])
{
    struct event e;
    // create mpi struct type for event
    const int event_property_count = 7;
    const int blocklengths[] = {4, 1, 1, 1, 1, 1, 1, 1};
    const MPI_Aint displacement[] = {0, 16, 24, 32, 40, 44, 48, 52};
    // printf("&e=%p, &e.num=%p, sizeof(e)=%lu\n", &e, &e.num, sizeof(e));
    const MPI_Datatype types[] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Datatype my_mpi_event_type;

    
    struct event event_recv_buff[X_SIZE * Y_SIZE+1];
    // printf();
	const int succeed_signal = 0;
	const int INTERNODE_COMM_TAG=0;
	const int BASE_COMM_TAG=1;
	const int SIMULATION_COMPLETED_SIGNAL = 2;
	int simulation_completion[X_SIZE * Y_SIZE+1] = {-1};
	MPI_Request base_comm_reqs[(X_SIZE * Y_SIZE+1)*2];
	// MPI_Request simulation_complete_reqs[X_SIZE * Y_SIZE+1];
    double tick, encryption_time, decryption_time, t0;
    int global_size, global_rank, size, rank, nneighbor;
    int random_num;
    // uint8_t* message;
    int*message;
    MPI_Comm node_comm;
    int r0, r1, recv_node;
    const int baserank_list[1] = {BASERANK};
    const int dim[2] = {X_SIZE, Y_SIZE}, period[2] = {0}, reorder = 0;
	int coord[2];
    int* neighbor_rank;
    MPI_Group global_group, node_group;
    MPI_Request *reqs;
    MPI_Request trash_req;

    MPI_Init(&argc, &argv);
    MPI_Type_create_struct(event_property_count, blocklengths, displacement, types, &my_mpi_event_type);
    MPI_Type_commit(&my_mpi_event_type);
    
    tick = MPI_Wtick();
    int* neighbor_result;
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
        srandom(time(NULL)+global_rank);

        message = calloc(MESSAGE_LEN, sizeof(int));
        MPI_Alloc_mem((MPI_Aint)(nneighbor * MESSAGE_LEN * sizeof(int)), MPI_INFO_NULL, &neighbor_result);
        MPI_Alloc_mem((MPI_Aint)(sizeof(MPI_Request) * 4), MPI_INFO_NULL, &reqs);
        MPI_Alloc_mem(nneighbor*sizeof(int), MPI_INFO_NULL, &neighbor_rank);

        // same as power(2, nbit), but more efficient
        int upperbound = 1 << N_BIT_RAND;
        int *num_recv = calloc(upperbound, sizeof(int));
        // e.occur_on_ranks = malloc(nneighbor);

        for (int current_i = 0; current_i < N_ITERATION; current_i++)
        {

            random_num = (int)(random() & ((1 << N_BIT_RAND) - 1));
            // random_num = 2;
            // printf("Random number is %d on rank %d coordinates are %d %d. Has %d neighbor\n", random_num, rank, coord[0], coord[1], nneighbor);

            // add padding of zeros following random num
            message[0] = random_num;
            // memcpy(message, &random_num, sizeof(int));

            t0 = MPI_Wtime(); 
            xor_encrypt(message, message, MESSAGE_LEN, key, KEY_SIZE);
            encryption_time = MPI_Wtime() - t0;

            for (int i=0, dim=0; dim<2; ++dim) {
                MPI_Cart_shift(node_comm, dim, 1, &r0, &r1);
                if (r0 >= 0) {
                    MPI_Isend(message, MESSAGE_LEN, MPI_INT, r0, INTERNODE_COMM_TAG, node_comm, &trash_req);
                    MPI_Irecv(&neighbor_result[i * MESSAGE_LEN], MESSAGE_LEN, MPI_INT, r0, INTERNODE_COMM_TAG, node_comm, &reqs[i]);
                    neighbor_rank[i] = r0;
                    i++;
                }
                if (r1 >= 0) {
                    MPI_Isend(message, MESSAGE_LEN, MPI_INT, r1, INTERNODE_COMM_TAG, node_comm, &trash_req);
                    MPI_Irecv(&neighbor_result[i * MESSAGE_LEN], MESSAGE_LEN, MPI_INT, r1, INTERNODE_COMM_TAG, node_comm, &reqs[i]);
                    neighbor_rank[i] = r1;
                    i++;
                }
            }
            MPI_Waitall(nneighbor, reqs, MPI_STATUSES_IGNORE);


            t0 = MPI_Wtime();
            xor_decrypt(neighbor_result, neighbor_result, nneighbor * MESSAGE_LEN, key, KEY_SIZE);
            decryption_time = MPI_Wtime() - t0;
			// for (int i=0;i<nneighbor;i++) { 
			// 	printf("rank %d received %d from rank %d\n", rank, (int)neighbor_result[i*MESSAGE_LEN], neighbor_rank[i]);
			// }

            for (int i=0;i <nneighbor; i++) {
                num_recv[(int)neighbor_result[i*MESSAGE_LEN]] ++; 
            }
            for (int i=0;i < upperbound; i++) {
				// printf("rank %d, num_recv[%d]=%d\n", rank, i, num_recv[i]);
                if (num_recv[i] >= RAND_TH)
                {
                    e.num = i;
                    e.n_times = num_recv[i];
                    e.iteration = current_i;
                    e.reference_rank = rank;
                    e.encryption_time = encryption_time;
                    e.decryption_time = decryption_time;
                    e.timestamp = time(NULL);
                    for (int j = 0, k=0; j < nneighbor; j++)
                    {
                        if (neighbor_result[j * MESSAGE_LEN] == i)
                        {
                            e.occur_on_ranks[k++] = neighbor_rank[j];
                        }
                    }
                    // print_event(e);
                    // print_event(e);
                    // printf("%d sent by %d\n", e.num, rank);
                    
                    MPI_Send(&e, 1, my_mpi_event_type, BASERANK, BASE_COMM_TAG, MPI_COMM_WORLD);
                    // print_event(e);
                    printf("event num %d being sent by rank %d\n", e.num, global_rank);
                    for (int j=i;j<upperbound;j++) {
                        num_recv[j] = 0;
                    }
                    break;
                }
                // clean up for next iteration
                num_recv[i] = 0;
            }
            usleep(1000 * INTERVAL);
        }
		MPI_Send(&succeed_signal, 1, MPI_INT, BASERANK, SIMULATION_COMPLETED_SIGNAL, MPI_COMM_WORLD);
		MPI_Comm_free(&node_comm);
        free(message);
    } else {
        // struct event p;
        // MPI_Recv(&p, 1, my_mpi_event_type, 2, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // // print_event(p);
        // printf("test elem %p\n", &p.reference_rank);
        // int* t = &p.num;
        // MPI_Recv(&p, 1, my_mpi_event_type, 5, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // printf("test elem %p\n", &p.num);
        // printf("t=%d\n",*t);
        // print_event(p);
        for (int i=0;i<global_size;i++) {
			if (i!=BASERANK) {
				MPI_Irecv(&event_recv_buff[i], 1, my_mpi_event_type, i, BASE_COMM_TAG, MPI_COMM_WORLD, &base_comm_reqs[i]);
				MPI_Irecv(&simulation_completion[i], 1, MPI_INT, i, SIMULATION_COMPLETED_SIGNAL, MPI_COMM_WORLD, &base_comm_reqs[i+global_size]);
			} else {
				base_comm_reqs[i] = MPI_REQUEST_NULL;
				base_comm_reqs[i+global_size] = MPI_REQUEST_NULL;
			}
		}
		int current_recv_rank;
        int tmpflag;
		int simulation_all_completed_flag = 0;
		MPI_Status status[X_SIZE*Y_SIZE+1];
        MPI_Status single_status;
        int remain_rank_num = global_size-1;

        while (!simulation_all_completed_flag) {
            MPI_Waitany(global_size*2, base_comm_reqs, &current_recv_rank, &single_status);
            if (current_recv_rank >= global_size) {
                printf("rank %d completed\n", current_recv_rank-global_size);
                base_comm_reqs[current_recv_rank-global_size] = MPI_REQUEST_NULL;
                base_comm_reqs[current_recv_rank] = MPI_REQUEST_NULL;
                remain_rank_num -= 1;
                if (remain_rank_num==0) {
                    simulation_all_completed_flag = 1;
                }
            } else {
                // printf("rank %d\n", current_recv_rank);
                // print_event(event_recv_buff[current_recv_rank]);
                printf("event num %d detected on rank %d\n", event_recv_buff[current_recv_rank].num, current_recv_rank);
                MPI_Irecv(&event_recv_buff[current_recv_rank], 1, my_mpi_event_type, current_recv_rank, BASE_COMM_TAG, MPI_COMM_WORLD, &base_comm_reqs[current_recv_rank]);
            }
		}
    }
    MPI_Finalize();
    return(0);
}


