
all: run
# X_SIZE and Y_SIZE control number of nodes in x dimension and y dimension
X_SIZE?=4
Y_SIZE?=5

# N_BIT_RAND defines how many bits random number generated occpy, controls the probability of event will decrease as it increase. The maximum value is 8 as the type of random number is uint8_t
N_BIT_RAND?=2

# N_INTERATION defines how many iterations the simulation will run
N_ITERATION?=200

# Each sensor will wail <INTERVAL> millionseconds between iterations
INTERVAL?=1

# MESSAGE_LEN represent the length of message that is passed among sensors. In plaintext, its first byte contains random number generated, and folloing bytes are 0 padding. It would be best to keep this as multiple of 16 (block size in CTR mode)
MESSAGE_LEN?=1024

# turn on debugging will result in information being printed to stdout
DEBUG?=0


wsn: wsn.c WjCryptLib_Aes.c WjCryptLib_AesCtr.c
	mpicc -DX_SIZE=$(X_SIZE) -DY_SIZE=$(Y_SIZE) -DDEBUG=$(DEBUG) -DN_BIT_RAND=$(N_BIT_RAND) -DINTERVAL=$(INTERVAL) -DN_ITERATION=$(N_ITERATION) -DMESSAGE_LEN=$(MESSAGE_LEN) -o $@ $^ -I. -lm -fopenmp
run: wsn
	mpirun -n "$$(( $(X_SIZE) * $(Y_SIZE) + 1 ))" --oversubscribe ./wsn
debug: wsn
	mpirun -n "$$(( $(X_SIZE) * $(Y_SIZE) + 1 ))"  xterm -e lldb -o run ./wsn.c  ./wsn 
clean:
	rm -rf wsn *.mpiP  logfile*.txt


