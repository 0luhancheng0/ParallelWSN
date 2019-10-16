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
#include <omp.h>
#include "WjCryptLib_AesCtr.h"

#ifndef DEBUG
#define DEBUG 0
#endif

struct event
{
    int occur_on_ranks[4];
//	int coordinate[2];
    double encryption_time;
    double decryption_time;
    long int timestamp;
    int n_times;
    int iteration;
    int reference_rank;
    uint8_t num;
};

struct node_summary
{
    int coordinate[2];
    int global_rank;
	int local_rank;	
	int threads_num;
	
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
	
	const uint8_t key[16] = {(uint8_t)0x16, (uint8_t)0x12, (uint8_t)0x75, (uint8_t)0x3e, (uint8_t)0xa8, (uint8_t)0xb7, (uint8_t)0x84, (uint8_t)0xe2, (uint8_t)0x1a, (uint8_t)0x2b, (uint8_t)0x3c, (uint8_t)0xdb, (uint8_t)0x9d, (uint8_t)0xcf, (uint8_t)0x4f, (uint8_t)0x3c};
	const uint8_t IV[8] = {(uint8_t)0x2b, (uint8_t)0x7e, (uint8_t)0x15, (uint8_t)0x16, (uint8_t)0x28, (uint8_t)0xae, (uint8_t)0xd2, (uint8_t)0xa6};
	AesCtrContext ctx;
	AesCtrInitialiseWithKey(&ctx, key, 16, IV);

    struct event e;
	// printf("sizeof e is %lu\n", sizeof(e));
	
	struct node_summary sensor_summary;
	MPI_Datatype mpi_sensor_summary_type;
	

    struct event *all_events;
    MPI_Status single_status;
    struct event event_recv_buff[X_SIZE * Y_SIZE+1];
	const int succeed_signal = 0;
	const int INTERNODE_COMM_TAG=0;
	const int BASE_COMM_TAG=1;
	const int SIMULATION_COMPLETED_SIGNAL = 2;
	int simulation_completion[X_SIZE * Y_SIZE+1] = {-1};
	MPI_Request base_comm_reqs[(X_SIZE * Y_SIZE+1)*2];
    double tick, encryption_time, decryption_time, t0;
    int global_size, global_rank, size, rank, nneighbor;
    uint8_t random_num;
    uint8_t* message;
    MPI_Comm node_comm;
    int r0, r1, recv_node;
    const int baserank_list[1] = {BASERANK};
    const int dim[2] = {X_SIZE, Y_SIZE}, period[2] = {0}, reorder = 0;
	int coord[2];
    int neighbor_rank[4];
	int upperbound = 1 << N_BIT_RAND;
    MPI_Group global_group, node_group;
    MPI_Request *reqs;
    MPI_Request trash_req;

    MPI_Init(&argc, &argv);
    MPI_Type_contiguous(5, MPI_INT, &mpi_sensor_summary_type);
    MPI_Type_commit(&mpi_sensor_summary_type);

    // MPI_Type_create_struct(event_property_count, blocklengths, displacement, types, &my_mpi_event_type);
    // MPI_Type_commit(&my_mpi_event_type);
 
   
    tick = MPI_Wtick();
    uint8_t* neighbor_result;
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
	#if DEBUG
	printf("i'm rank %d\n", global_rank);
	#endif
    
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

		// create sensor summary
		sensor_summary.global_rank = global_rank;
        sensor_summary.local_rank = rank;
        memcpy(&sensor_summary.coordinate[0], &coord[0], sizeof(int) * 2);
        sensor_summary.threads_num = omp_get_num_threads();
        


        // setup different random number seed on different rank
        srandom(time(NULL)+global_rank);

		// memory allocation
        message = calloc(MESSAGE_LEN, 1);
		#if DEBUG
		printf("rank %d has %d neighbor\n", global_rank, nneighbor);
		#endif
        MPI_Alloc_mem((MPI_Aint)(nneighbor * MESSAGE_LEN), MPI_INFO_NULL, &neighbor_result);
        MPI_Alloc_mem((MPI_Aint)(sizeof(MPI_Request) * 4), MPI_INFO_NULL, &reqs);
        // MPI_Alloc_mem(nneighbor * sizeof(int), MPI_INFO_NULL, &neighbor_rank);


        // same as power(2, nbit), but more efficient
        int *num_recv = calloc(upperbound, sizeof(int));

        for (int current_i = 0; current_i < N_ITERATION; current_i++)
        {
            random_num = (uint8_t)(random() & ((1 << N_BIT_RAND) - 1));
			#if DEBUG
			random_num=0;
			printf("random number is %d rank %d\n", random_num, global_rank);
			#endif 

            // add padding of zeros following random num
            message[0] = random_num;

			// encryption
            t0 = MPI_Wtime(); 
            // xor_encrypt(message, message, MESSAGE_LEN * sizeof(int));
			AesCtrXor(&ctx, message, message, MESSAGE_LEN);
			AesCtrSetStreamIndex(&ctx, 0);
            encryption_time = MPI_Wtime() - t0;

			// inter node communication
            for (int i=0, dim=0; dim<2; ++dim) {
                MPI_Cart_shift(node_comm, dim, 1, &r0, &r1);
                if (r0 >= 0) {
                    MPI_Isend(message, MESSAGE_LEN, MPI_UINT8_T, r0, INTERNODE_COMM_TAG, node_comm, &trash_req);
                    MPI_Irecv(&neighbor_result[i * MESSAGE_LEN], MESSAGE_LEN, MPI_UINT8_T, r0, INTERNODE_COMM_TAG, node_comm, &reqs[i]);
                    neighbor_rank[i] = r0;
                    i++;
                }
                if (r1 >= 0) {
                    MPI_Isend(message, MESSAGE_LEN, MPI_UINT8_T, r1, INTERNODE_COMM_TAG, node_comm, &trash_req);
                    MPI_Irecv(&neighbor_result[i * MESSAGE_LEN], MESSAGE_LEN, MPI_UINT8_T, r1, INTERNODE_COMM_TAG, node_comm, &reqs[i]);
                    neighbor_rank[i] = r1;
                    i++;
                }
            }
            MPI_Waitall(nneighbor, reqs, MPI_STATUSES_IGNORE);
			#if DEBUG
			printf("rank %d finished internode comm in iteration %d\n", global_rank, current_i);
			#endif

			// decryption
            t0 = MPI_Wtime();

			// count number of occurence for each number
			// AesCtrXor(&ctx, neighbor_result, neighbor_result, MESSAGE_LEN*nneighbor);
            for (int i=0;i <nneighbor; i++) {
				// AesCtrXor(&ctx, neighbor_result, neighbor_result, MESSAGE_LEN);
				AesCtrXor(&ctx, &neighbor_result[i*MESSAGE_LEN], &neighbor_result[i*MESSAGE_LEN], MESSAGE_LEN);
				AesCtrSetStreamIndex(&ctx, 0);
                num_recv[(int)neighbor_result[i*MESSAGE_LEN]] ++; 
				#if DEBUG	
				printf("current pointer %llu on rank %d\n", ctx.StreamIndex, global_rank);
				printf("on rank %d, num_recv[%d]=%d\n",global_rank, i, num_recv[(int)neighbor_result[i*MESSAGE_LEN]]);
				#endif
            }
            decryption_time = MPI_Wtime() - t0;
			
			// send event to base station
            for (int i=0;i < upperbound; i++) {
                if (num_recv[i] >= RAND_TH)
                {
                    e.num = i;
					// e.coordinate[0] = coord[0];
					// e.coordinate[1] = coord[1];
                    e.n_times = num_recv[i];
                    e.iteration = current_i;
                    e.reference_rank = global_rank;
                    e.encryption_time = encryption_time * (double)num_recv[i];
                    e.decryption_time = decryption_time;
                    e.timestamp = time(NULL);
                    for (int j = 0, k=0; j < nneighbor; j++)
                    {
                        if (neighbor_result[j * MESSAGE_LEN] == i)
                        {
                            e.occur_on_ranks[k++] = neighbor_rank[j];
                        }
                    }
                    AesCtrXor(&ctx, &e, &e,sizeof(e));
                    AesCtrSetStreamIndex(&ctx, 0);
					#if DEBUG
					struct event *p = malloc(sizeof(e));
					AesCtrXor(&ctx, &e, p,sizeof(e));
					AesCtrSetStreamIndex(&ctx, 0);
					//print_event(*p);

					#endif
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
        AesCtrXor(&ctx, &sensor_summary, &sensor_summary, sizeof(int) * 5);
        // AesCtrSetStreamIndex(&)
        MPI_Send(&sensor_summary, 1, mpi_sensor_summary_type, BASERANK, SIMULATION_COMPLETED_SIGNAL, MPI_COMM_WORLD);
        printf("rank %d completed\n", global_rank);
		// MPI_Send(&succeed_signal, 1, MPI_UINT8_T, BASERANK, SIMULATION_COMPLETED_SIGNAL, MPI_COMM_WORLD);
    } 
		else {

        int max_possible_events = X_SIZE*Y_SIZE*N_ITERATION;
        int event_storage_p=0;
        int total_event_num = 0;
        struct node_summary* all_sensor_summary;
        all_sensor_summary = malloc(sizeof(struct node_summary) * global_size);
        all_events = malloc(max_possible_events * sizeof(struct event));

		// spawn set of receive for nodes
        for (int i = 0; i < global_size; i++)
        {
            if (i!=BASERANK) {

				// this receive is for actual event message
                MPI_Irecv(&event_recv_buff[i], sizeof(struct event), MPI_UINT8_T, i, BASE_COMM_TAG, MPI_COMM_WORLD, &base_comm_reqs[i]);
                // this receive is for the signal at the end of communication
                MPI_Irecv(&all_sensor_summary[i], 1, mpi_sensor_summary_type, i, SIMULATION_COMPLETED_SIGNAL, MPI_COMM_WORLD, &base_comm_reqs[i + global_size]);
            } else {
				// set request as MPI_REQUEST_NULL for root node itself
				base_comm_reqs[i] = MPI_REQUEST_NULL;
				base_comm_reqs[i+global_size] = MPI_REQUEST_NULL;
			}
		}
		int current_recv_rank;
		int simulation_all_completed_flag = 0;
        int remain_rank_num = global_size-1;
        double encryption_time = 0;
        double decryption_time = 0;
        int *event_activation_num = calloc(upperbound, sizeof(int));

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
                AesCtrXor(&ctx, &all_sensor_summary[current_recv_rank-global_size], &all_sensor_summary[current_recv_rank-global_size], sizeof(int) * 5);
                AesCtrSetStreamIndex(&ctx, 0);
                printf("global %d, local %d coord (%d %d)\n", all_sensor_summary[current_recv_rank - global_size].local_rank, all_sensor_summary[current_recv_rank - global_size].global_rank, all_sensor_summary[current_recv_rank - global_size].coordinate[0], all_sensor_summary[current_recv_rank - global_size].coordinate[1]);
                // all_sensor_summary[current_recv_rank-global_size].global_rank = current_recv_rank-global_size;

                // if event is received
            } else {
				// store received event
				AesCtrXor(&ctx, &event_recv_buff[current_recv_rank], &all_events[event_storage_p], sizeof(e));
				AesCtrSetStreamIndex(&ctx, 0);
                
                total_event_num++;
                encryption_time += all_events[event_storage_p].encryption_time;
                decryption_time += all_events[event_storage_p].decryption_time;
				event_activation_num[(int)all_events[event_storage_p].num] ++;
				event_storage_p++;

				// respawn receive
                MPI_Irecv(&event_recv_buff[current_recv_rank], sizeof(struct event), MPI_UINT8_T, current_recv_rank, BASE_COMM_TAG, MPI_COMM_WORLD, &base_comm_reqs[current_recv_rank]);
            }
		}

		// create logfile
		time_t now = time(NULL);
		struct tm* localt = localtime(&now);
		
		char filename[50];
		int tmp = sprintf(filename, "logfile_");
		tmp += strftime(&filename[tmp], 50-tmp, "%H_%M_%S", localt);
        // filename = asctime(localtime(&now));
		sprintf(&filename[tmp], ".txt");
        FILE *f = fopen(filename, "w");

		// the header of logfile contains an overview of network
        char header[800] = "Event detection in a fully distributed wireless sensor network - WSN\n\n";
        int header_p = strlen(header);
        header_p += sprintf(header + header_p, "network configuration overview: \n");
        header_p += sprintf(header + header_p, "Simulation will run %d iterations with %d milliseconds time interval between each pair of consecutive iteration\n", N_ITERATION, INTERVAL);
        header_p += sprintf(header + header_p, "Network have %d nodes in X dimension and %d nodes in Y dimension\n", X_SIZE, Y_SIZE);
        header_p += sprintf(header + header_p, "The size of random number is %d bits, which indicate event number is bound by range [0, %d)\n", N_BIT_RAND, 1 << N_BIT_RAND);
        header_p += sprintf(header + header_p, "process and topology summary: \n");
        for (int i=0; i<global_size;i++) {
            if (i != BASERANK) {
                printf("%d %d %d %d %d\n",all_sensor_summary[i].global_rank, all_sensor_summary[i].local_rank, all_sensor_summary[i].coordinate[0], all_sensor_summary[i].coordinate[1], all_sensor_summary[i].threads_num);
                header_p += sprintf(header + header_p, "\tprocess rank %d in MPI_COMM_WORLD has local rank %d in sensors communicatior. The coordinate is (%d, %d).%d OpenMP threads are available on this process\n", all_sensor_summary[i].global_rank, all_sensor_summary[i].local_rank, all_sensor_summary[i].coordinate[0], all_sensor_summary[i].coordinate[1], all_sensor_summary[i].threads_num);
            }
        }
        
        header_p += sprintf(header + header_p, "details of communication:\n");
        header_p += sprintf(header + header_p, "event activation summary: \n");
		for (int i=0;i<upperbound;i++) {
			header_p += sprintf(header + header_p, "\tevent %d is activated %d times\n", i, event_activation_num[i]);
		}


        header_p += sprintf(header + header_p, "For each activation: encryption time : %lf seconds\ndecryption time : %lf seconds\n", encryption_time, decryption_time);
        header_p += sprintf(header + header_p, "number of message pass between base station and nodes : %d\n", total_event_num + X_SIZE * Y_SIZE);
        header_p += sprintf(header + header_p, "number of message passing happened among nodes: %d\n", (X_SIZE * (Y_SIZE - 1) + Y_SIZE * (X_SIZE - 1)) * 2 * N_ITERATION);
        header_p += sprintf(header + header_p, "total events detected: %d\n\n", total_event_num);
        fwrite(header, sizeof(char), header_p, f);
        free(&header);

        // write log for received event details
        char *event_details_log = malloc(total_event_num*sizeof(char)*200);
        int event_p = 0;
        int perv_iteration = -1;
        for (int i=0;i<total_event_num;i++) {
            if (all_events[i].iteration != perv_iteration) {
                event_p += sprintf(event_details_log + event_p, "Iteration %d\n\n", all_events[i].iteration);
                perv_iteration = all_events[i].iteration;
            }
            event_p += sprintf(event_details_log + event_p, "event number %d detected on local rank %d (rank %d globally) with coordinate (%d, %d)\n", all_events[i].num, all_sensor_summary[all_events[i].reference_rank].local_rank, all_sensor_summary[all_events[i].reference_rank].global_rank, all_sensor_summary[all_events[i].reference_rank].coordinate[0], all_sensor_summary[all_events[i].reference_rank].coordinate[1]);

            event_p += sprintf(event_details_log + event_p, "Timestamp : %s", asctime(localtime(&all_events[i].timestamp)));
            event_p += sprintf(event_details_log + event_p, "adjacent nodes are (local): ");
            for (int j=0;j<all_events[i].n_times;j++) {
                event_p += sprintf(event_details_log + event_p, "%d ", all_events[i].occur_on_ranks[j]);
            }
            event_p += sprintf(event_details_log + event_p, "\n");

            // add more stuff here
            event_p += sprintf(event_details_log + event_p, "\n");
        }

        
        fwrite(event_details_log, sizeof(char), event_p, f);
        
        fclose(f);


    }
	#if DEBUG
	printf("rank %d is going to end\n", global_rank);
	#endif 
    MPI_Finalize();
    return(0);
}


