// length of key in bytes
#define KEY_SIZE 16

// 128 bit key
uint8_t key[16] =  { (uint8_t) 0x2b, (uint8_t) 0x7e, (uint8_t) 0x15, (uint8_t) 0x16, (uint8_t) 0x28, (uint8_t) 0xae, (uint8_t) 0xd2, (uint8_t) 0xa6, (uint8_t) 0xab, (uint8_t) 0xf7, (uint8_t) 0x15, (uint8_t) 0x88, (uint8_t) 0x09, (uint8_t) 0xcf, (uint8_t) 0x4f, (uint8_t) 0x3c };

void encrypt(uint8_t* plaintext, uint8_t* ciphertext, int message_length, uint8_t* k, int k_size);
void decrypt(uint8_t* ciphertext, uint8_t* plaintext, int message_length, uint8_t* k, int k_size);
