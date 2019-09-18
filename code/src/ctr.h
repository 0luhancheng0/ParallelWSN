
#define KYESIZE 128
#define BLOCKLEN 4
#include <stdint.h>
typedef struct
{
    uint8_t RoundKey[KYESIZE];
    uint8_t Iv[BLOCKLEN];
} ctx;
void encryption(struct AES_ctx *ctx, uint8_t *buf, uint32_t length)