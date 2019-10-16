// size of each dimension

#ifndef X_SIZE
#define X_SIZE 4
#endif
#ifndef Y_SIZE
#define Y_SIZE 5
#endif
// the rank of base station in MPI_COMM_WORLD
#ifndef BASERANK
#define BASERANK 0
#endif
// the number of bits each random number have, which also defines the upper bound of random number as 2^n-1, maximum 8 bits
#ifndef N_BIT_RAND
#define N_BIT_RAND 2
#endif
// run N_INTERATION number times
#ifndef N_ITERATION
#define N_ITERATION 200
#endif
// run sleep ITERVAL millisecond each iteration
#ifndef INTERVAL
#define INTERVAL 1
#endif
// length of message
#ifndef MESSAGE_LEN
#define MESSAGE_LEN 1024
#endif
// threshold for random number, event will be reported when occurence is greater or equal to this number
#ifndef RAND_TH
#define RAND_TH 1
#endif




