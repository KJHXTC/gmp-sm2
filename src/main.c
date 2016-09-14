/*
 * test.c
 *
 *  Created on: 2012-3-19
 *      Author: nacle
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "gmp.h"
#include "ec_operations.h"
#include "sm2_func.h"
#include "debug.h"
#include "yl_base_tools.h"
//#define   str_P   "21A6F0638B79118292EFF1194966C79EBEEC66CEF9DE0F16179B404141E8CE264E0C0C4E00D3"
//#define   str_P   "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF"
#define str_P  "8542D69E4C044F18E8B92435BF6FF7DE457283915C45517D722EDB8B08F1DFC3"
//#define str_A 	 "21A6F0638B79118292EFF1194966C79EBEEC66CEF9DE0F16179B404141E8CE264E0C0C4E00D0"
//#define   str_A   "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC"
#define str_A  "787968B4FA32C3FD2417842E73BBFEFF2F3C848B6831D7E0EC65228B3937E498"
//#define   str_B  "17CE4670BA97C022AF1EBB13DDD273E477550721413D3B9AB11C2CF998FFA43CA34286D4AD26"
//#define   str_B	"28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93"
#define str_B  "63E4C6D3B23B0C849CF84241484BFE48F61D59A5B16BA06E6E12D1DA27C5249A"
//#define   str_N	"FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123"
#define str_N  "8542D69E4C044F18E8B92435BF6FF7DD297720630485628D5AE74EE7C32E79B7"
//#define   str_GX  "06128FE16E46D271E62228E0D7910089970D206CB36AD33AC83228B43222B610489F51C0EAAD"
//#define   str_GX  "32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7"
#define str_GX "421DEBD61B62EAB6746434EBC3CC315E32220B3BADD50BDC4C4E6C147FEDD43D"
//#define   str_GY  "0DD684C6FC67575354153DEB60E37DA054878F0175CE90A0EC8560FF1C7D225CBFF2FA763A7C"
//#define   str_GY  "BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0"
#define str_GY "0680512BCBB42C07D47349D2153B70C4E5D7FDFCBFA36EA1A85841B9E46E09A2"
//priveta key
//#define   str_D   "F06AC3DC3ED41ACA03EED04CDB3CAF38ED54766E2A38CFE6E1BF0CD04845B434"
//#define str_D  "128B2FA8BD433C6C068C8D803DFF79792A519A55171B1B650C23661D15897263"
#define str_D  "1649AB77A00637BD5E2EFE283FBF353534AA7F7CB89463F208DDBC2920BB0DA0"
//randmo
#define   str_R   "34914C20251A59A2C311102944C600430A02285A0433144228142A1848004C14"
//#define str_R  "6CB28D99385C175C94F94E934817663FC176D925DD72B727260DBAAE1FB2F96F"
//hase value
//#define   str_E  "01EAA84A92E56558E7336AF1F93695EFE5E50A8F66B4B6C710D99A2DD93C910E52119F"
//#define   str_E   "23f8ad84dcc3277eb4050b8ee1f5965b"
//#define str_E  "B524F552CD82B8B028476E005C377FB19A87E6FC682D48BB5D42E3D9B9EFFE76"
#define str_E 	"656E6372797074696F6E207374616E64617264"
//Pubkey
//#define str_Pkx	 	"8940911ba05e78858f209dfa753dd4d3ad14c76d071315c9d890e373f61cf12a"
//#define str_Pky  	"8ce9232ea9d449537ec29a2032f257077018bfe14abe3bb887e40a9ab76874a5"
//#define str_Pkx	 "0AE4C7798AA0F119471BEE11825BE46202BB79E2A5844495E97C04FF4DF2548A"
//#define str_Pky  "7C0240F88F1CD4E16352A73C17B7F16F07353E53A176D684A9FE0C6BB798E857"
#define str_Pkx	 	"435B39CCA8F3B508C1488AFC67BE491A0F7BA07E581A0E4849A5CF70628A7E0A"
#define str_Pky  	"75DDBA78F15FEECB4C7895E2C1CDF5FE01DEBB2CDBADF45399CCF77BBA076A42"
//signed value
#define  str_RSM  "c22fffa3c2a421ae387dea9515a0bd9d3238b8df376964baa6f47f0de76e6420"
#define  str_SSM  "b3e037ecdff14f4facf99149dacbeb2feff129ede29e242231c7df9e40872d23"
//#define  str_RSM  "40F1EC59F793D9F49E09DCEF49130D4194F79FB1EED2CAA55BACDB49C4E755D1"
//#define  str_SSM  "6FC6DAC32C5D5CF10C77DFB20F7C2EB667A457872FB09EC56327A67EC7DEEBE7"

gmp_randstate_t stat;
int main(){
	int ret ,i;
	int base =16;
	mpz_t curv[2];
	mpz_t base_point[2];
	mpz_t Q[2];
	mpz_t result[2];
	mpz_t e;
	mpz_t p;
	mpz_t n;
	mpz_t r;
	mpz_t s;
	mpz_t d;
	mpz_t seed;
	unsigned char * e_char ,*orgtest;
	unsigned long reslen;
	unsigned char  e_str[100];
    int sd , clen;
    int klen = 152;
    clen = VBYTE*3+klen/8+1;
    unsigned char C[clen];
	mpz_init(seed);
	mpz_init(result[0]);
	mpz_init(result[1]);
	mpz_init_set_str(curv[0],str_A ,base);
	mpz_init_set_str(curv[1],str_B ,base);
	mpz_init_set_str(base_point[0],str_GX ,base);
	mpz_init_set_str(base_point[1],str_GY ,base);
	mpz_init_set_str(Q[0],str_Pkx ,base);
	mpz_init_set_str(Q[1],str_Pky ,base);
	mpz_init_set_str(e,str_E ,base);
	mpz_init_set_str(p,str_P ,base);
	mpz_init_set_str(n,str_N ,base);
	mpz_init_set_str(r,str_RSM ,base);
	mpz_init_set_str(s,str_SSM ,base);
	mpz_init_set_str(d,str_D ,base);
#ifdef DEBUG_PRINT
	gmp_printf("curv[0] is %Zx\n" , curv[0]);
	gmp_printf("curv[1] is %Zx\n" , curv[1]);
	gmp_printf("base_point[0] is %Zx\n" , base_point[0]);
	gmp_printf("base_point[1] is %Zx\n" , base_point[1]);
	gmp_printf("Q[0] is %Zx\n" , Q[0]);
	gmp_printf("Q[1] is %Zx\n" , Q[1]);
	gmp_printf("e is %Zx\n" , e);
	gmp_printf("p is %Zx\n" , p);
	gmp_printf("n is %Zx\n" , n);
	gmp_printf("r is %Zx\n" , r);
	gmp_printf("s is %Zx\n" , s);
#endif
	gmp_randinit(stat, GMP_RAND_ALG_LC, 120);
	srand( (unsigned) getpid());
	sd=rand();
	mpz_set_ui(seed, sd);
	gmp_randseed(stat, seed);
	mpz_get_str(e_str,16,e);
	e_char = hextochs(e_str);
	pubkey_encryption(curv,base_point,&p,&n,Q,e_char,klen,C,&reslen);
#ifdef DEBUG_PRINT
	printf("C = ");
	for(i=0 ;i<reslen;i++){
		printf("%02x",*(C+i));
	}
	printf("\n");
#endif
	orgtest = (unsigned char*)malloc(sizeof(unsigned char)*reslen);
	ret = privkey_decryption(curv,base_point,&p,&n,&d,C,klen,orgtest);
	if(ret== -1)
		printf("privkey_decryption fail!\n");
#ifdef DEBUG_PRINT
	printf("orgtest = ");
	for(i=0 ;i<klen/8;i++){
		printf("%02x",*(orgtest+i));
	}
	printf("\n");
#endif
}
