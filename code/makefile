
all: run
# X and Y control number of nodes in x dimension and y dimension
X?=4
Y?=5

# max number of OpenMP threads
OMP_NUM_THREADS?=1

# N_BIT_RAND defines how many bits random number generated occpy, controls the probability of event will decrease as it increase. The maximum value is 8 as the type of random number is uint8_t
N_BIT_RAND?=3

# N_INTERATION defines how many iterations the simulation will run
N_ITERATION?=200

# Each sensor will wail <INTERVAL> milliseconds between iterations
INTERVAL?=10

# MESSAGE_LEN represent the length of message that is passed among sensors. In plaintext, its first byte contains random number generated, and folloing bytes are 0 padding. It would be best to keep this as multiple of 16 (block size in CTR mode)
MESSAGE_LEN?=1024

# turn on debugging will result in information being printed to stdout
DEBUG?=0

# print the message before and after encryption 
PRINT_MESSAGE?=0

# extra flags passed to mpi run
FLAGS?=

wsn: wsn.c WjCryptLib_Aes.c WjCryptLib_AesCtr.c
	mpicc -DX_SIZE=$(X) -DY_SIZE=$(Y) -DDEBUG=$(DEBUG) -DN_BIT_RAND=$(N_BIT_RAND) -DINTERVAL=$(INTERVAL) -DN_ITERATION=$(N_ITERATION) -DMESSAGE_LEN=$(MESSAGE_LEN) -DPRINT_MESSAGE=$(PRINT_MESSAGE)  -o $@ $^ -I. -lm -fopenmp
run: wsn
	mpirun -n "$$(( $(X) * $(Y) + 1 ))" -x OMP_NUM_THREADS=$(OMP_NUM_THREADS) $(FLAGS) ./wsn
debug: wsn
	mpirun -n "$$(( $(X) * $(Y) + 1 ))" --oversubscribe xterm -e lldb -o run ./wsn.c  ./wsn 
clean:
	rm -rf wsn *.mpiP  logfile*.txt a.out slurm*.out smux*.out


