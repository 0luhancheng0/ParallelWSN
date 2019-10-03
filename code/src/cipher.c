#include <stdint.h>
#include <stdio.h>

void xor_encrypt(uint8_t* plaintext, uint8_t* ciphertext, int message_length, const uint8_t* key, int keysize) {
    int chunk_size = 1, nthreads, thread_id, i, j;
    #pragma omp parallel shared(nthreads, chunk_size, ciphertext) private(thread_id, i, j)
    {
        #pragma omp for schedule(static, chunk_size)
        for (i=0;i<message_length;i++)
        {
            ciphertext[i] = plaintext[i] ^ key[i % keysize];
        }
        
    }
 
};
void xor_decrypt(uint8_t* ciphertext, uint8_t* plaintext, int message_length, const uint8_t* key, int keysize) {
    xor_encrypt(ciphertext, plaintext, message_length, key, keysize);
};

