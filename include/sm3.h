/*
 * sm3.h
 *
 * 为使此算法兼容32位、64位下Linux或Windows系统，
 * 选择 int 来表示 32 位整数。
 * 消息长度最大限定为 2**32 - 1（单位：比特），
 * 且为 8 的倍数（消息的最小单元为字节）。
 */
#ifndef _SM3_H_
#define _SM3_H_

/*
 * SM3算法产生的哈希值大小（单位：字节）
 */
#define SM3_HASH_SIZE 32 

/*
 * SM3上下文
 */
typedef struct SM3Context
{
    unsigned int intermediateHash[SM3_HASH_SIZE / 4];
    unsigned char messageBlock[64];
} SM3Context;

/*
 * SM3计算函数
 */
unsigned char *SM3Calc(const unsigned char *message, 
        unsigned int messageLen, unsigned char digest[SM3_HASH_SIZE]);

/**
 * \brief          Output = SM3( input buffer )
 *
 * \param input    buffer holding the  data
 * \param ilen     length of the input data
 * \param output   SM3 checksum result
 */
void sm3( unsigned char *input, int ilen,
           unsigned char output[32]);

#endif // _SM3_H_