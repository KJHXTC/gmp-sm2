/*
 * base_tools.h
 */


#ifndef YL_BASE_TOOLS_H_
#define YL_BASE_TOOLS_H_
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#define VBITS 256
#define VBYTE ((VBITS+7)/8)
unsigned char *hextochs ( const char* ascii );
unsigned char *chstohex ( const unsigned char * chs );
void  sm2Kdf(unsigned char *z,unsigned int zlen,unsigned int klen,unsigned char *output);
#endif /* YL_BASE_TOOLS_H_ */
