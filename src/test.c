/*
 * test.c
 *
 *  Created on: 2012-3-19
 *      Author: nacle
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gmp.h"
#include "sm2_func.h"
#include "yl_base_tools.h"


gmp_randstate_t stat;
int main(){
	
	
	
	for (int i=1; i<100; i++)
	{
		test(i);
	}
	
	
	return 0;
}

void test(int cnt){
	printf("%d *************************************************\n\n", cnt);
	int i;
	char* pkx= "CBC6D43900D8257C9D643792FBD9C597BFD676A9ED1BE7668A14CD9A8CFB86E9";
	char* pky= "36572DBB0C3996933F8F4677158A6C07F5B3BCEFE472016ECE04599E37240C86";
	char* data = "656e6372797074696f6e207374616e64617264";
	UBYTE* result=NULL;
	ULONG msglen = strlen(data)/2;
	ULONG len = sm2_encrypt(256, pkx, pky, data, msglen, &result);
	
	printf("test.result %x= [", result);
	for(i=0; i<len; i++)
	{
		printf("%02X", result[i]);
	}	
	printf("]\n");
	printf("free result::%X\n", (result));
	free(result);
	result=NULL;
	printf("test done\n\n\n\n\n");
}