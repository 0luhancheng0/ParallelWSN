#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "test.h"
#include "cipher.h"
#include <memory.h>
int main(int argc, char* argv[]) {
    srand(time(NULL));
    uint8_t random_num = (uint8_t)(rand() & ((1 << N_BIT_RAND) - 1));
    printf("size of randn %d\n", sizeof(random_num));
    uint8_t* message;
    message = (uint8_t*)calloc(MESSAGE_LEN, sizeof(uint8_t));
    memcpy(message, &random_num, sizeof(random_num));
    printf("message length is %d\n", MESSAGE_LEN);
    printf("message is %d\n", message[1]);
    free(message);
    return(0);
}

void padding(uint8_t randn, uint8_t* out) {

}
