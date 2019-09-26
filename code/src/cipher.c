#include <stdint.h>
#include <omp.h>
#include <stdio.h>

void encrypt(uint8_t* plaintext, uint8_t* ciphertext, int message_length, uint8_t* key, int keysize) {
    int chunk_size = 1, nthreads, thread_id, i, j;
    #pragma omp parallel shared(nthreads, chunk_size, ciphertext) private(thread_id, i, j)
    {
        nthreads = omp_get_num_threads();
        thread_id = omp_get_thread_num();
        #pragma omp for schedule(static, chunk_size)
        for (i=0;i<message_length;i++)
        {
            ciphertext[i] = plaintext[i] ^ key[i % keysize];
        }
        
    }
 
};
void decrypt(uint8_t* ciphertext, uint8_t* plaintext, int message_length, uint8_t* key, int keysize) {
    encrypt(ciphertext, plaintext, message_length, key, keysize);
};

