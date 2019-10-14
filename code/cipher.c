#include <stdint.h>
#include <stdio.h>
// length of key in bytes
#define KEY_SIZE 16

// 128 bit key
static const uint8_t key[16] = {(uint8_t)0x2b, (uint8_t)0x7e, (uint8_t)0x15, (uint8_t)0x16, (uint8_t)0x28, (uint8_t)0xae, (uint8_t)0xd2, (uint8_t)0xa6, (uint8_t)0xab, (uint8_t)0xf7, (uint8_t)0x15, (uint8_t)0x88, (uint8_t)0x09, (uint8_t)0xcf, (uint8_t)0x4f, (uint8_t)0x3c};
// static const int key[4] = {123,432,5435, 43453};
// IV

void xor_encrypt(uint8_t* plaintext, uint8_t* ciphertext, int message_length) {
	// plaintext = (uint8_t*)plaintext;
	// ciphertext = (uint8_t*)ciphertext;
	
    int chunk_size = 1, nthreads, thread_id, i, j;
    #pragma omp parallel shared(nthreads, chunk_size, ciphertext) private(thread_id, i, j)
    {
        #pragma omp for schedule(static, chunk_size)
        for (i=0;i<message_length;i++)
        {
            ciphertext[i] = plaintext[i] ^ key[i % KEY_SIZE];
        }
        
    }
 
};
void xor_decrypt(uint8_t* ciphertext, uint8_t* plaintext, int message_length) {
    xor_encrypt(ciphertext, plaintext, message_length);
};

