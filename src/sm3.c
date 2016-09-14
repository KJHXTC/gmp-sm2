/*
 * sm3.c
 */
#include <stdio.h>
#include <memory.h>
#include "sm3.h"

/*
 * 判断运行环境是否为小端
 */
static const int endianTest = 1;
#define IsLittleEndian() (*(char *)&endianTest == 1)

/*
 * 向左循环移位
 */
#define LeftRotate(word, bits) ( (word) << (bits) | (word) >> (32 - (bits)) )

/*
 * 反转四字节整型字节序
 */
unsigned int *ReverseWord(unsigned int *word)
{
    unsigned char *byte, temp;

    byte = (unsigned char *)word;
    temp = byte[0];
    byte[0] = byte[3];    
    byte[3] = temp;

    temp = byte[1];
    byte[1] = byte[2];
    byte[2] = temp;
    return word;
}

/*
 * T
 */
unsigned int T(int i)
{
    if (i >= 0 && i <= 15)
        return 0x79CC4519;
    else if (i >= 16 && i <= 63)
        return 0x7A879D8A;
    else
        return 0;
}

/*
 * FF
 */
unsigned int FF(unsigned int X, unsigned int Y, unsigned int Z, int i)
{
    if (i >= 0 && i <= 15)
        return X ^ Y ^ Z;
    else if (i >= 16 && i <= 63)
        return (X & Y) | (X & Z) | (Y & Z);
    else
        return 0;
}

/*
 * GG
 */
unsigned int GG(unsigned int X, unsigned int Y, unsigned int Z, int i)
{
    if (i >= 0 && i <= 15)
        return X ^ Y ^ Z;
    else if (i >= 16 && i <= 63)
        return (X & Y) | (~X & Z);
    else
        return 0;
}

/*
 * P0
 */
unsigned int P0(unsigned int X)
{
    return X ^ LeftRotate(X, 9) ^ LeftRotate(X, 17);
}

/*
 * P1
 */
unsigned int P1(unsigned int X)
{
    return X ^ LeftRotate(X, 15) ^ LeftRotate(X, 23);
}

/*
 * 初始化函数
 */
void SM3Init(SM3Context *context)
{
    context->intermediateHash[0] = 0x7380166F;
    context->intermediateHash[1] = 0x4914B2B9;
    context->intermediateHash[2] = 0x172442D7;
    context->intermediateHash[3] = 0xDA8A0600;
    context->intermediateHash[4] = 0xA96F30BC;
    context->intermediateHash[5] = 0x163138AA;
    context->intermediateHash[6] = 0xE38DEE4D;
    context->intermediateHash[7] = 0xB0FB0E4E;
}

/*
 * 处理消息块
 */
void SM3ProcessMessageBlock(SM3Context *context)
{
    int i;
    unsigned int W[68];
    unsigned int W_[64];
    unsigned int A, B, C, D, E, F, G, H, SS1, SS2, TT1, TT2;

    /* 消息扩展 */
    for (i = 0; i < 16; i++)
    {
        W[i] = *(unsigned int *)(context->messageBlock + i * 4);
        if (IsLittleEndian())
            ReverseWord(W + i);
        //printf("%d: %x\n", i, W[i]);    
    }
    for (i = 16; i < 68; i++)
    {
        W[i] = P1(W[i - 16] ^ W[i - 9] ^ LeftRotate(W[i - 3], 15)) 
                ^ LeftRotate(W[i - 13], 7) 
                ^ W[i - 6];
        //printf("%d: %x\n", i, W[i]);    
    }
    for (i = 0; i < 64; i++)
    {
        W_[i] = W[i] ^ W[i + 4];
        //printf("%d: %x\n", i, W_[i]);    
    }

    /* 消息压缩 */
    A = context->intermediateHash[0];
    B = context->intermediateHash[1];
    C = context->intermediateHash[2];
    D = context->intermediateHash[3];
    E = context->intermediateHash[4];
    F = context->intermediateHash[5];
    G = context->intermediateHash[6];
    H = context->intermediateHash[7];
    for (i = 0; i < 64; i++)
    {
        SS1 = LeftRotate((LeftRotate(A, 12) + E + LeftRotate(T(i), i)), 7);
        SS2 = SS1 ^ LeftRotate(A, 12);
        TT1 = FF(A, B, C, i) + D + SS2 + W_[i];
        TT2 = GG(E, F, G, i) + H + SS1 + W[i];
        D = C;
        C = LeftRotate(B, 9);
        B = A;
        A = TT1;
        H = G;
        G = LeftRotate(F, 19);
        F = E;
        E = P0(TT2);
    }
    context->intermediateHash[0] ^= A;
    context->intermediateHash[1] ^= B;
    context->intermediateHash[2] ^= C;
    context->intermediateHash[3] ^= D;
    context->intermediateHash[4] ^= E;
    context->intermediateHash[5] ^= F;
    context->intermediateHash[6] ^= G;
    context->intermediateHash[7] ^= H;
}

/*
 * SM3算法主函数
 */
unsigned char *SM3Calc(const unsigned char *message, 
        unsigned int messageLen, unsigned char digest[SM3_HASH_SIZE])
{
    SM3Context context;
    unsigned int i, remainder, bitLen;

    /* 初始化上下文 */
    SM3Init(&context);

    /* 对前面的消息分组进行处理 */
    for (i = 0; i < messageLen / 64; i++)
    {
        memcpy(context.messageBlock, message + i * 64, 64);
        SM3ProcessMessageBlock(&context);
    }

    /* 填充消息分组，并处理 */
    bitLen = messageLen * 8;
    if (IsLittleEndian())
        ReverseWord(&bitLen);
    remainder = messageLen % 64;
    memcpy(context.messageBlock, message + i * 64, remainder);
    context.messageBlock[remainder] = 0x80;
    if (remainder <= 55)
    {
        /* 长度按照大端法占8个字节，该程序只考虑长度在 2**32 - 1（单位：比特）以内的情况，
         * 故将高 4 个字节赋为 0 。*/
        memset(context.messageBlock + remainder + 1, 0, 64 - remainder - 1 - 8 + 4);
        memcpy(context.messageBlock + 64 - 4, &bitLen, 4);    
        SM3ProcessMessageBlock(&context);
    }
    else
    {
        memset(context.messageBlock + remainder + 1, 0, 64 - remainder - 1);
        SM3ProcessMessageBlock(&context);
        /* 长度按照大端法占8个字节，该程序只考虑长度在 2**32 - 1（单位：比特）以内的情况，
         * 故将高 4 个字节赋为 0 。*/
        memset(context.messageBlock, 0, 64 - 4);
        memcpy(context.messageBlock + 64 - 4, &bitLen, 4);    
        SM3ProcessMessageBlock(&context);
    }

    /* 返回结果 */
    if (IsLittleEndian())
        for (i = 0; i < 8; i++)
            ReverseWord(context.intermediateHash + i);
    memcpy(digest, context.intermediateHash, SM3_HASH_SIZE);

    return digest;
}

/*
 * output = SM3( input buffer )
 */
void sm3( unsigned char *input, int ilen,
           unsigned char output[32] )
{
    SM3Context ctx;

    SM3Calc(input, ilen, output);
}
