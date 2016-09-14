/*
 * sm2_func.h
 *
 *  Created on: 2012-3-19
 *      Author: nacle
 */
#ifndef SM2_FUNC
#define SM2_FUNC

#define UBYTE unsigned char
#define ULONG unsigned long

// SM2 国密参数
#define	str_P	"FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF"
#define	str_A   "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC"
#define	str_B   "28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93"
#define	str_N   "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123"
#define	str_GX  "32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7" 
#define	str_GY  "BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0"

extern gmp_randstate_t stat;
/* 返回加密数据结果一般 1+96+消息长度就是结果长度 */
ULONG sm2_encrypt(
		int keylen, 
		const char* pkx,		// 公钥 X 16进制字符串表示
		const char* pky,		// 公钥 Y 16进制字符串表示
		const char* message,	// 待加密消息 16进制字符串表示
		const ULONG datalen,	// 消息长度，转换为实际字节的长度，一般是 16进制串长度的一半
		UBYTE ** result			// 加密结果 PC||C1||C3||C2
		);

void pubkey_encryption(
		mpz_t * curv,
		mpz_t * base_point,
		mpz_t * p,
		mpz_t * n,
		mpz_t * Pb,
		UBYTE * msg, ULONG klen,
		UBYTE * res, ULONG *reslen
		);

#endif
