#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <memory.h>
#include "test.h"
#include "cipher.h"

int main(int argc, char* argv[]) {
    srand(time(NULL));
    uint8_t random_num = (uint8_t)(rand() & ((1 << N_BIT_RAND) - 1));
    uint8_t* message;
    uint8_t* cipher;
    message = (uint8_t*)calloc(MESSAGE_LEN, sizeof(uint8_t));
    cipher = (uint8_t*)malloc(MESSAGE_LEN * sizeof(uint8_t));
    memcpy(message, &random_num, sizeof(random_num));
    encrypt(message, cipher, MESSAGE_LEN, key, KEY_SIZE);
    decrypt(cipher, cipher, MESSAGE_LEN, key, KEY_SIZE);
    int res = memcmp(message, cipher, MESSAGE_LEN);
    free(message);
    return(res);
}


