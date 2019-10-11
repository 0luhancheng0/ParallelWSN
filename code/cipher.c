#include <stdint.h>
#include <stdio.h>

void xor_encrypt(int* plaintext, int* ciphertext, int message_length, const int* key, const int keysize) {
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
void xor_decrypt(int* ciphertext, int* plaintext, int message_length, const int* key, const int keysize) {
    xor_encrypt(ciphertext, plaintext, message_length, key, keysize);
};

