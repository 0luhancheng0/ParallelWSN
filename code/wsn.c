#include <stdio.h>
#include <memory.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include <unistd.h>
#include "wsn.h"
#include "cipher.h"
#include <omp.h>

struct event
{
    int occur_on_ranks[4];
	int coordinate[2];
    double encryption_time;
    double decryption_time;
    long int timestamp;
    int n_times;
    int iteration;
    int num;
    int reference_rank;
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
    struct event *all_events;
    MPI_Status single_status;
	
    // const int event_property_count = 10;
    // const int blocklengths[] = {4, 2, 1, 1, 1, 1, 1, 1, 1, 1};
    // const MPI_Aint displacement[] = {0, 16, 24, 32, 40, 48, 52, 56, 60, sizeof(e)};
    // const MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_UB};

    // MPI_Datatype my_mpi_event_type;
   
    struct event event_recv_buff[X_SIZE * Y_SIZE+1];
	const int succeed_signal = 0;
	const int INTERNODE_COMM_TAG=0;
	const int BASE_COMM_TAG=1;
	const int SIMULATION_COMPLETED_SIGNAL = 2;
	int simulation_completion[X_SIZE * Y_SIZE+1] = {-1};
	MPI_Request base_comm_reqs[(X_SIZE * Y_SIZE+1)*2];
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
 
    // MPI_Type_create_struct(event_property_count, blocklengths, displacement, types, &my_mpi_event_type);
    // MPI_Type_commit(&my_mpi_event_type);
 
   
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

		// get coordinate in cart
        MPI_Cart_coords(node_comm, rank, 2, coord);
        get_neighbor_count(coord, &nneighbor);

        // setup different random number seed on different rank
        srandom(time(NULL)+global_rank);

		// memory allocation
        message = calloc(MESSAGE_LEN, sizeof(int));
        MPI_Alloc_mem((MPI_Aint)(nneighbor * MESSAGE_LEN * sizeof(int)), MPI_INFO_NULL, &neighbor_result);
        MPI_Alloc_mem((MPI_Aint)(sizeof(MPI_Request) * 4), MPI_INFO_NULL, &reqs);
        MPI_Alloc_mem(nneighbor*sizeof(int), MPI_INFO_NULL, &neighbor_rank);

        // same as power(2, nbit), but more efficient
        int upperbound = 1 << N_BIT_RAND;
        int *num_recv = calloc(upperbound, sizeof(int));

        for (int current_i = 0; current_i < N_ITERATION; current_i++)
        {
            random_num = (int)(random() & ((1 << N_BIT_RAND) - 1));

            // add padding of zeros following random num
            message[0] = random_num;

			// encryption
            t0 = MPI_Wtime(); 
            xor_encrypt(message, message, MESSAGE_LEN * sizeof(int));
            encryption_time = MPI_Wtime() - t0;

			// inter node communication
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

			// decryption
            t0 = MPI_Wtime();
            xor_decrypt(neighbor_result, neighbor_result, nneighbor * MESSAGE_LEN * sizeof(int));
            decryption_time = MPI_Wtime() - t0;

			// count number of occurence for each number
            for (int i=0;i <nneighbor; i++) {
                num_recv[(int)neighbor_result[i*MESSAGE_LEN]] ++; 
            }
			
			// send event to base station
            for (int i=0;i < upperbound; i++) {
                if (num_recv[i] >= RAND_TH)
                {
                    e.num = i;
					e.coordinate[0] = coord[0];
					e.coordinate[1] = coord[1];
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
                    xor_encrypt(&e, &e, sizeof(e));
                    MPI_Send(&e, sizeof(e), MPI_UINT8_T, BASERANK, BASE_COMM_TAG, MPI_COMM_WORLD);
					// clean up for next iteration
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

		// send succeed signal to base station when simulation complete
		MPI_Send(&succeed_signal, 1, MPI_INT, BASERANK, SIMULATION_COMPLETED_SIGNAL, MPI_COMM_WORLD);
		MPI_Comm_free(&node_comm);
        free(message);
    } else {
        int max_possible_events = X_SIZE*Y_SIZE*N_ITERATION;
        int event_storage_p=0;
        int total_event_num = 0;
        all_events = malloc(max_possible_events * sizeof(struct event));

		// spawn set of receive for nodes
        for (int i = 0; i < global_size; i++)
        {
            if (i!=BASERANK) {

				// this receive is for actual event message
				// MPI_Irecv(&event_recv_buff[i], 1, my_mpi_event_type, i, BASE_COMM_TAG, MPI_COMM_WORLD, &base_comm_reqs[i]);
                MPI_Irecv(&event_recv_buff[i], sizeof(struct event), MPI_UINT8_T, i, BASE_COMM_TAG, MPI_COMM_WORLD, &base_comm_reqs[i]);
                // this receive is for the signal at the end of communication
                MPI_Irecv(&simulation_completion[i], 1, MPI_INT, i, SIMULATION_COMPLETED_SIGNAL, MPI_COMM_WORLD, &base_comm_reqs[i + global_size]);
            } else {
				// set request as MPI_REQUEST_NULL for root node itself
				base_comm_reqs[i] = MPI_REQUEST_NULL;
				base_comm_reqs[i+global_size] = MPI_REQUEST_NULL;
			}
		}
		int current_recv_rank;
		int simulation_all_completed_flag = 0;
        int remain_rank_num = global_size-1;
        double average_encryption_time = 0;
        double average_decryption_time = 0;

        while (!simulation_all_completed_flag) {

			// both event receiver and signal receiver are stored in the same array, so they can be waited at the same line
            MPI_Waitany(global_size*2, base_comm_reqs, &current_recv_rank, &single_status);

			// if complete signal is received
            if (current_recv_rank >= global_size) {

				// set both receiver to null to disable further waiting
                base_comm_reqs[current_recv_rank-global_size] = MPI_REQUEST_NULL;
                base_comm_reqs[current_recv_rank] = MPI_REQUEST_NULL;
                remain_rank_num -= 1;

				// flag if all nodes are finished
                if (remain_rank_num==0) {
                    simulation_all_completed_flag = 1;
                }
			
			// if event is received
            } else {
				// store received event
                xor_decrypt(&event_recv_buff[current_recv_rank], &event_recv_buff[current_recv_rank], sizeof(e));
                memcpy(&all_events[event_storage_p++], &event_recv_buff[current_recv_rank], sizeof(struct event));
                
                total_event_num++;
                average_encryption_time += event_recv_buff[current_recv_rank].encryption_time;
                average_decryption_time += event_recv_buff[current_recv_rank].decryption_time / (double)event_recv_buff[current_recv_rank].n_times;

				// respawn receive
                MPI_Irecv(&event_recv_buff[current_recv_rank], sizeof(struct event), MPI_UINT8_T, current_recv_rank, BASE_COMM_TAG, MPI_COMM_WORLD, &base_comm_reqs[current_recv_rank]);
            }
		}
        average_encryption_time = average_encryption_time / total_event_num;
        average_decryption_time = average_decryption_time / total_event_num;

		// create logfile
        char filename[] = "log.txt";
        FILE *f = fopen(filename, "w");

		// the header of logfile contains an overview of network
        char header[500] = "Event detection in a fully distributed wireless sensor network - WSN\n\n";
        ;
        int header_p = strlen(header);
        header_p += sprintf(header + header_p, "network configuration overview: \n");
        header_p += sprintf(header + header_p, "Simulation will run %d iterations with %d milliseconds time interval between each pair of consecutive iteration\n", N_ITERATION, INTERVAL);
        header_p += sprintf(header + header_p, "Network have %d nodes in X dimension and %d nodes in Y dimension\n", X_SIZE, Y_SIZE);
        header_p += sprintf(header + header_p, "The size of random number is %d bits, which indicate event number is bound by range [0, %d)\n", N_BIT_RAND, 1 << N_BIT_RAND);
        header_p += sprintf(header + header_p, "For each message: average encryption time : %lf seconds\naverage decryption time : %lf seconds\n", average_encryption_time, average_decryption_time);
        header_p += sprintf(header + header_p, "number of message pass between base station and nodes : %d\n", total_event_num + X_SIZE * Y_SIZE);
        header_p += sprintf(header + header_p, "number of message passing happened among nodes: %d\n", (X_SIZE * (Y_SIZE - 1) + Y_SIZE * (X_SIZE - 1)) * 2 * N_ITERATION);
        header_p += sprintf(header + header_p, "total events detected: %d\n\n", total_event_num);
        header_p += sprintf(header + header_p, "details of nodes involved in communication:\n");

		// write log for received event details
        char *event_details_log = malloc(total_event_num*sizeof(char)*300);
        int event_p = 0;
        int perv_iteration = 0;
        for (int i=0;i<total_event_num;i++) {
            if (all_events[i].iteration != perv_iteration) {
                event_p += sprintf(event_details_log + event_p, "Iteration %d\n\n", all_events[i].iteration);
                perv_iteration = all_events[i].iteration;
            }
            event_p += sprintf(event_details_log + event_p, "event number %d detected on rank (local) %d with coordinate (%d, %d)\n", all_events[i].num, all_events[i].reference_rank, all_events[i].coordinate[0], all_events[i].coordinate[1]);
            event_p += sprintf(event_details_log + event_p, "Timestamp : %s", asctime(localtime(&all_events[i].timestamp)));
            event_p += sprintf(event_details_log + event_p, "adjacent nodes are : ");
            for (int j=0;j<all_events[i].n_times;j++) {
                event_p += sprintf(event_details_log + event_p, "%d ", all_events[i].occur_on_ranks[j]);
            }
            event_p += sprintf(event_details_log + event_p, "\n");

            // add more stuff here
            event_p += sprintf(event_details_log + event_p, "\n");
        }

        fwrite(header, sizeof(char), header_p, f);
        fwrite(event_details_log, sizeof(char), event_p, f);
        
        fclose(f);


    }
    MPI_Finalize();
    return(0);
}


