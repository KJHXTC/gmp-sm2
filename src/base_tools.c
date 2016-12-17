/*
 * base_tools.c
 */
#include "base_tools.h"
#include "sm3.h"
void sm2Kdf(unsigned char *z,unsigned int zlen,unsigned int klen,unsigned char *output)
{
	unsigned int hlen, outlen,i;
	unsigned char *ha,*Zn;
	unsigned char hai[VBYTE];
	
	hlen = (klen+ VBITS-1)/VBITS;

	ha=(unsigned char *)malloc(VBYTE*(hlen+1));
	Zn=(unsigned char *)malloc(zlen+4);
	memset(ha, 0, VBYTE*hlen);
	memset(Zn, 0, zlen+4);
	memcpy(Zn,z,zlen);
	for(i=1;i<= hlen;i++)
	{
		char tmp[8+1]={0};
		sprintf(tmp, "%08X", i);
		unsigned char* byte4 ;
		byte4 = hextochs(tmp);
		
		memset(Zn+zlen+0, *(byte4+0), 1);
		memset(Zn+zlen+1, *(byte4+1), 1);
		memset(Zn+zlen+2, *(byte4+2), 1);
		memset(Zn+zlen+3, *(byte4+3), 1);

		free(byte4);
		sm3(Zn,zlen+4,hai);
		memcpy(ha+(i-1)*VBYTE, hai,VBYTE);  
	}
	
	if ((klen % VBITS) > 0)
	{
		memcpy(ha+(i-1)*VBYTE, hai, VBYTE);
	}

	outlen=(klen+7)/8;
	memcpy(output, ha, outlen);
	
	free(ha);
	free(Zn);
}



unsigned char *hextochs (const char* ascii )
{
    int len = strlen(ascii) ;
    if( len%2 != 0 )
    {
		char hexs[99999+1]={0};
		hexs[0] = '0';
		if (strlen(ascii)<99999)
		{
			return hextochs(strcat(hexs, ascii));
		}else{
			return hextochs(strncat(hexs, ascii, 99999));
		}
        
    }
    unsigned char *chs = NULL ;
    chs = (unsigned char*)calloc( len / 2 + 1, sizeof(unsigned char) );                // calloc chs
    int  i = 0 ;
    char ch[2] = {0};
    while( i < len )
    {
        ch[0] = ( (int)ascii[i] > 64 ) ? ( ascii[i]%16 + 9 ) : ascii[i]%16 ;
        ch[1] = ( (int)ascii[i + 1] > 64 ) ? ( ascii[i + 1]%16 + 9 ) : ascii[i + 1]%16 ;

        chs[i/2] = (unsigned char)( ch[0]*16 + ch[1] );
        i += 2;
    }

    return chs ;            // chs ����ǰδ�ͷ�
}


unsigned char  *chstohex (const unsigned char * chs )
{
	unsigned char  hex[16] = { '0', '1', '2', '3', '4', '5', '6', \
        '7', '8','9', 'A', 'B', 'C', 'D', 'E', 'F'};

    int len = strlen ( chs );
    unsigned char * ascii = NULL ;
    ascii = (unsigned char *)calloc ( len * 2 + 1, sizeof(unsigned char ) );            // calloc ascii
    int i = 0;
    while( i < len )
    	{
        ascii[i*2] = hex[(int)( (unsigned char )chs[i] / 16 )] ;
        ascii[i*2 + 1] = hex[(int)( (unsigned char )chs[i] % 16 )] ;
        ++i;
    	}

    return ascii;                    // ascii ����֮ǰδ�ͷ�
}

