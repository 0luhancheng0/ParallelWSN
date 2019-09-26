// the rank of base station in MPI_COMM_WORLD
#define BASERANK 0

// the number of bits each random number have, which also defines the upper bound of random number as 2^n-1, maximum 1byte
#define N_BIT_RAND 4

// run N_INTERATION number times
#define N_ITERATION 10000

// run sleep ITERVAL millisecond each iteration
#define INTERVAL 100

// length of message in bytes after padding precending 0s must be greater than N_BIT_RAND / 8, need to be multiple of keysize
#define MESSAGE_LEN 512
