
#include <stdint.h>
// xor_encrypt and xor_decrypt are the same function, which xor key with first argument then store result in second argument
void xor_encrypt(void* plaintext, void* ciphertext, int message_length);
void xor_decrypt(void* ciphertext, void* plaintext, int message_length);
