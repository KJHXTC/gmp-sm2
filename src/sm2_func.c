/*
 * sm2_func.c
 *
 *  Created on: 2012-3-19
 *      Author: nacle
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "debug.h"
#include "gmp.h"
#include "int_arithmetic.h"
#include "ec_operations.h"
#include "sm2_func.h"
#include "sm3.h"
#include "yl_base_tools.h"

/*
 *use pubkey pb encrypt the message m, the result is res .
 * */
void pubkey_encryption(
		mpz_t * curv, mpz_t * base_point,mpz_t * p, mpz_t * n, 
		mpz_t * Pb,
		unsigned char *m, unsigned long mlen,
		unsigned char * res ,unsigned long *reslen){

	gmp_randinit_default(stat);
		
	mpz_t r,k,k1;
	mpz_t p1[2];
	mpz_t p2[2];

	mpz_init(r);
	mpz_init(k);
	mpz_init(k1);
	mpz_init(p1[0]);
	mpz_init(p1[1]);
	mpz_init(p2[0]);
	mpz_init(p2[1]);
	
	int	i,kchar,len;

	kchar = mlen;
	unsigned char K[kchar];
	unsigned char c1[VBYTE*2];
	unsigned char c2[kchar];
	unsigned char c3[VBYTE];
	unsigned char tmp[VBYTE*2];
	unsigned char *xmy;
	unsigned char *px;
	unsigned char *py;

	printf("1. 产生随机数k \n");
//1. 产生随机数k k in range(1, N)
	while (mpz_sgn(r) == 0) {
		myprimegenerator(&k1);
		printf("1.1 myprimegenerator \n");
		mpz_mod(k, k1, *n);
		gmp_printf("k1 is %ZX\n", k1);

		while (mpz_sgn(k) == 0){
			myprimegenerator(&k1);
		}
			
		printf("2.计算K倍的G 点 \n");
		point_mult_window(curv, base_point, &k, p1, p);
		mpz_mod(r, p1[0], *n);
	}
	gmp_printf("x1 is %ZX\n" ,p1[0]);
	gmp_printf("y1 is %ZX\n" ,p1[1]);



	// 计算得到C1
	
	px = hextochs(mpz_get_str(NULL, 16, p1[0]));
	printf("colloc c1 px: ptr::%x\n", px);
	py = hextochs(mpz_get_str(NULL, 16, p1[1]));
	printf("colloc c1 py: ptr::%x\n", py);
	printf("3.计算C1\n");
	memcpy(c1, px, VBYTE);
	memcpy(c1+VBYTE, py, VBYTE);
	printf("free c1 px: ptr::%x\n", px);
	free(px);
	printf("free c1 py: ptr::%x\n", py);
	free(py);
#ifdef DEBUG_PRINT
	printf("C1 %d = [", VBYTE);
	for(i=0 ;i<(VBYTE*2);i++){
                printf("%02X",*(c1+i));
        }
	printf("]\n");
#endif
	// 公钥计算是否无穷点S
	printf("4.计算Pb的k倍点C2 \n");
	point_mult_window(curv, Pb, &k, p2, p);

	px = hextochs(mpz_get_str(NULL, 16, p2[0]));
	printf("colloc c2 py: ptr::%x\n", px);
	py = hextochs(mpz_get_str(NULL, 16, p2[1]));
	printf("colloc c2 py: ptr::%x\n", py);
	memcpy(tmp, px, VBYTE);
	memcpy(tmp+VBYTE, py, VBYTE);
	printf("5.产生密钥数据K \n");
	sm2Kdf(tmp, VBYTE*2, mlen*8, K);
	printf("t %d = [", kchar);
	for(i=0 ;i<kchar;i++){
		printf("%02X", K[i]);
	}
	printf("]\n");
	printf("m %d = [", kchar);
	for(i=0 ;i<kchar;i++){
		printf("%02X", m[i]);
	}
	printf("]\n");
	for (i = 0; i < kchar; i++)
		c2[i] = m[i] ^ K[i];
#ifdef DEBUG_PRINT
	gmp_printf("x2 is %ZX\n" ,p2[0]);
	gmp_printf("y2 is %ZX\n" ,p2[1]);

	printf("C2 %d = [", kchar);
	for(i=0 ;i<kchar;i++){
		printf("%02X",*(c2+i));
	}
	printf("]\n");
#endif
//	计算C3
	xmy=(unsigned char *)malloc(VBYTE*2+kchar);
	memset(xmy,0, (VBYTE*2)+kchar);
	memcpy(xmy, px, (VBYTE));
	memcpy(xmy+ VBYTE, m, kchar);
	memcpy(xmy+ VBYTE+kchar, py, (VBYTE));
	len = VBYTE*2 + kchar;
	sm3(xmy,len,c3);

#ifdef DEBUG_PRINT
	printf("C3 %d = [", VBYTE);
	for(i=0 ;i<VBYTE;i++){
		printf("%02X",*(c3+i));
	}
	printf("]\n");
#endif
// 将C1 C3 C2 打入res
	res[0]=(unsigned char) 0x04;
	int offset=1;
	memcpy(res +offset, c1, VBYTE*2);
	offset += VBYTE*2;
	memcpy(res +offset, c3, VBYTE);
	offset += VBYTE;
	memcpy(res +offset, c2, mlen);
	offset += mlen;
	*reslen = offset;
	
#ifdef DEBUG_PRINT
	printf("result %d = [", offset);
	for(i=0 ;i<offset;i++){
		printf("%02X",*(res+i));
	}
	printf("]\n");
#endif
	printf("free c2 px: ptr::%x\n", px);
	//free(px);
	printf("free c2 py: ptr::%x\n", py);
	//free(py);
	printf("free: ptr::%x\n", xmy);
	free(xmy);


	mpz_clear(r);
	mpz_clear(k);
	mpz_clear(k1);

	gmp_randclear (stat);

	printf("clear r k k1 done;;;\n");
	mpz_clear(p1[0]);
	printf("clear p1x done;;;\n");
	mpz_clear(p1[1]);
	printf("clear p1y done;;;\n");
	mpz_clear(p2[0]);
	mpz_clear(p2[1]);
	printf("clear p2 done;;;\n");
	printf("clear done;;;\n");

}


ULONG sm2_encrypt(int keylen, const char* pkx, const char* pky, const char* message, const ULONG datalen, UBYTE ** result)
{
	
	int base =16;
	mpz_t curv[2];
	mpz_t base_point[2];
	mpz_t Q[2];
	mpz_t p;
	mpz_t n;

	
	// 设置 公共参数 
	mpz_init_set_str(curv[0], str_A, base);
	mpz_init_set_str(curv[1], str_B, base);
	mpz_init_set_str(base_point[0], str_GX, base);
	mpz_init_set_str(base_point[1], str_GY, base);
	mpz_init_set_str(p, str_P, base);
	mpz_init_set_str(n, str_N, base);

	// 设置公钥
	mpz_init_set_str(Q[0], pkx, base);
	mpz_init_set_str(Q[1], pky, base);
	//推测长度 flag +C1x+C1y+C3 +M 
	ULONG clen;
	clen = 1+ VBYTE*3+datalen;
	*result = (UBYTE*) malloc(clen);
/*
	gmp_printf("curv[0] is %Zx\n" , curv[0]);
	gmp_printf("curv[1] is %Zx\n" , curv[1]);
	gmp_printf("base_point[0] is %Zx\n" , base_point[0]);
	gmp_printf("base_point[1] is %Zx\n" , base_point[1]);
	gmp_printf("p is %Zx\n" , p);
	gmp_printf("n is %Zx\n" , n);
*/
#ifdef DEBUG_PRINT
	gmp_printf("Q[0] is %Zx\n" , Q[0]);
	gmp_printf("Q[1] is %Zx\n" , Q[1]);
	gmp_printf("message is %s\n" , message);
#endif

	ULONG reslen=0;

	UBYTE * data ;
	data = hextochs(message);
	gmp_printf("call pkenc data\n");
	pubkey_encryption(curv, base_point, &p, &n, Q, data, datalen, *result, &reslen);
	gmp_printf("end pkenc data\n");
//	释放指针
	free(data);

	mpz_clear(curv);
	mpz_clear(base_point);
	mpz_clear(p);
	mpz_clear(n);
	mpz_clear(Q);
	
	return reslen;
}