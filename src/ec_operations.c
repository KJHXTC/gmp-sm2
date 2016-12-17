// --------------------------------------------------------------------
//
//  File:        ec_operations.c
//  Date:        11/03
//  Last update: 04/10
//  Description: Basic operations on elliptic curves' group
//
//  (C) 2003, Elisavet Konstantinou & Yiannis Stamatiu & Christos Zaroliagis
//                 {konstane,stamatiu,zaro}@ceid.upatras.gr
//
//   Extended by: Bhanu Prakash & Pratik Poddar & Bernard Menezes 
//                 {prakashb,pratik, bernard}@cse.iitb.ac.in
// --------------------------------------------------------------------


# include <stddef.h>
# include <stdio.h>
# include <stdlib.h>
//# include <unistd.h>

# include "gmp.h"
# include "int_arithmetic.h"
# include "ec_operations.h"


/* generates randomly the a and b coefficients (curv[0] and curv[1]
   respectively) of the elliptic curve modulo the prime p */
void rand_curve(mpz_t * curv, mpz_t * p)
{
	mpz_t a1;
	mpz_t b1;

	mpz_init(a1);
	mpz_init(b1);

	myprimegenerator(&a1);
	myprimegenerator(&b1);

	mpz_mod(curv[0], a1, *p);
	mpz_mod(curv[1], b1, *p);

	mpz_clear(a1);
	mpz_clear(b1);
}


/* print the coefficients of an elliptic curve */
void print_curve(mpz_t * curv)
{
	printf(" Curve form: y^2=x^3+ax+b\n");

	printf("first curve parameter (a):");
	mpz_out_str(stdout, 10, curv[0]);
	printf("\n");

	printf("second curve parameter (b):");
	mpz_out_str(stdout, 10, curv[1]);
	printf("\n");

}

/* generate a random point base in the elliptic curve curv modulo p */
void rand_point(mpz_t *curv, mpz_t *p, mpz_t *base)
{
	mpz_t prime;
	mpz_t x1, h1, h2;

	mpz_init(prime);
	mpz_init(x1);
	mpz_init(h1);
	mpz_init(h2);

	mpz_urandomb(x1, stat, bitlength);

	do {
		mpz_mod(base[0], x1, *p);

		mpz_powm_ui(h1, base[0], 3, *p);
		myzmulmod(&h2, &curv[0], &base[0], p);
		myzaddmod(&h1, &h1, &h2, p);
		myzaddmod(&h1, &h1, &curv[1], p);

		if (mpz_sgn(h1) == 0)
			mpz_set_ui(base[1], 0);

		else
			myzsqrtmod(&base[1], &h1, p);


		mpz_add_ui(x1, x1, (long) 1);

	} while (mpz_sgn(base[1]) == 0);


	mpz_clear(prime);
	mpz_clear(x1);
	mpz_clear(h1);
	mpz_clear(h2);
}


/* print a point */
void print_point(mpz_t * point)
{
	printf("x coefficient:");
	mpz_out_str(stdout, 10, point[0]);
	printf("\n");

	printf("y coefficient:");
	mpz_out_str(stdout, 10, point[1]);
	printf("\n");

}


/* addition of two points p1 and p2 modulo p that belong in the
  elliptic curve curv. The result of the addition is the point p3 */
void add_point(mpz_t * curv, mpz_t * p1, mpz_t * p2, mpz_t * p3, mpz_t * p)
{
	mpz_t a1;
	mpz_t b1;
	mpz_t c1;
	mpz_t x1;

	mpz_t lamda;

	long flag = 0;

	mpz_init(a1);
	mpz_init(b1);
	mpz_init(c1);
	mpz_init(x1);

	mpz_init(lamda);

	// point at infinity = (0, 0)

	if ((mpz_sgn(p1[0]) == 0) && (mpz_sgn(p1[1]) == 0)) {
		mpz_set(p3[0], p2[0]);
		mpz_set(p3[1], p2[1]);
	}

	else if ((mpz_sgn(p2[0]) == 0) && (mpz_sgn(p2[1]) == 0)) {
		mpz_set(p3[0], p1[0]);
		mpz_set(p3[1], p1[1]);
	}

	else {
		if (mpz_cmp(p1[0], p2[0]) == (long) 0) {
			if (mpz_cmp(p1[1], p2[1]) == (long) 0)
			{
				if (mpz_sgn(p2[1]) == 0) {
					printf("point at infinity\n");
					mpz_set_ui(p3[0], 0);
					mpz_set_ui(p3[1], 0);
					flag = 1;
				}
				else
				{
					//x1 = x2 and p2!=-p1
					mpz_powm_ui(b1, p1[0], 2, *p);//b1 = x1^2modp
					myzsmulmod(&c1, &b1, (long) 3, p);//c1 = 3*x1^2 mod p
					myzaddmod(&b1, &c1, &curv[0], p);//b1 = c1+a mod p
					myzsmulmod(&a1, &p1[1], (long) 2, p);//a1 = 2*y1 mod p
					myzdivmod(&lamda, &b1, &a1, p);//lamda = b1/a1 mod p
				}
			}
			else
			{		// point at infinity
				mpz_set_ui(p3[0], 0);
				mpz_set_ui(p3[1], 0);
				flag = 1;
			}
		}
		else
		{
			//x1!=x2
			myzsubmod(&a1, &p2[1], &p1[1], p);//y2-y1
			myzsubmod(&b1, &p2[0], &p1[0], p);//x2-x1

			if ((mpz_sgn(a1) == (long) -1))
			{
				mpz_abs(a1, a1);
				if ((mpz_sgn(b1) == (long) -1)) {
					mpz_abs(b1, b1);
					myzdivmod(&lamda, &a1, &b1, p);//lamda = a1/b1 mod p
				} else {
					myzdivmod(&lamda, &a1, &b1, p);//lamda = a1/b1 mod p
					mpz_sub(lamda, *p, lamda);
				}
			}
			else if ((mpz_sgn(b1) == (long) -1))
			{
				mpz_abs(b1, b1);
				myzdivmod(&lamda, &a1, &b1, p);
				mpz_sub(lamda, *p, lamda);
			}

			else
				myzdivmod(&lamda, &a1, &b1, p);

		}

		mpz_powm_ui(a1, lamda, 2, *p);//a1 = lamda^2 mod p

		myzsubmod(&b1, &a1, &p1[0], p);//b1 = a1 - x1

		myzsubmod(&c1, &b1, &p2[0], p);//c1 = b1 -x2

		mpz_mod(x1, c1, *p);

		mpz_set(p3[0], x1);

		myzsubmod(&a1, &p1[0], &x1, p);//a1 = x1 -x3
		myzmulmod(&b1, &a1, &lamda, p);//b1 = lamda * a1
		myzsubmod(&c1, &b1, &p1[1], p);//c1 = b1-y1
		mpz_set(p3[1], c1);

	}
	if (flag == 1) {
		mpz_set_ui(p3[0], 0);
		mpz_set_ui(p3[1], 0);
	}


	mpz_clear(a1);
	mpz_clear(b1);
	mpz_clear(c1);
	mpz_clear(x1);

	mpz_clear(lamda);

}


/* implement scalar multiplication of the point p1 and k modulo p
   in the elliptic curve curv using the binary method */
void point_mult(mpz_t * curv, mpz_t * p1, mpz_t * k, mpz_t * result,
				mpz_t * p)
{
	mpz_t Y[2];
	mpz_t P[2];
	mpz_t R[2];
	mpz_t R1[2];

	long b = 0, i = 0;
	mpz_t k1;
	mpz_t k2;
	mpz_t b1;

	mpz_init(Y[0]);
	mpz_init(Y[1]);
	mpz_init(P[0]);
	mpz_init(P[1]);
	mpz_init(R[0]);
	mpz_init(R[1]);
	mpz_init(R1[0]);
	mpz_init(R1[1]);

	mpz_init(k1);
	mpz_init(k2);
	mpz_init(b1);


	mpz_set(P[0], p1[0]);
	mpz_set(P[1], p1[1]);


	mpz_set(k2, *k);

	while (mpz_sgn(k2) == 1) {

		i = i + 1;

		b = mpz_tdiv_qr_ui(k1, b1, k2, 2);

		if (b == (long) 1) {
			add_point(curv, Y, P, R1, p);
			mpz_set(Y[0], R1[0]);
			mpz_set(Y[1], R1[1]);

		}

		mpz_set(k2, k1);

		add_point(curv, P, P, R, p);
		mpz_set(P[0], R[0]);
		mpz_set(P[1], R[1]);

	}

	mpz_set(result[0], Y[0]);
	mpz_set(result[1], Y[1]);


	mpz_clear(Y[0]);
	mpz_clear(Y[1]);
	mpz_clear(P[0]);
	mpz_clear(P[1]);
	mpz_clear(R[0]);
	mpz_clear(R[1]);
	mpz_clear(R1[0]);
	mpz_clear(R1[1]);

	mpz_clear(k1);
	mpz_clear(k2);
	mpz_clear(b1);
}



/* implement scalar multiplication of the point p1 and k modulo p
   in the elliptic curve curv using the naf method */
void point_mult_naf(mpz_t * curv, mpz_t * p1, mpz_t * k, mpz_t * result,
				mpz_t * p)
{
	mpz_t Y[2];
	mpz_t P[2];
	mpz_t Q[2];
	mpz_t R[2];
	mpz_t R1[2];

	long b = 0, i = 0;
	mpz_t k1;
	mpz_t k2;
	mpz_t b1;
	mpz_t one;
	mpz_init_set_str(one,"1",10);
	mpz_init(Y[0]);
	mpz_init(Y[1]);
	mpz_init(P[0]);
	mpz_init(P[1]);
	mpz_init(Q[0]);
	mpz_init(Q[1]);
	mpz_init(R[0]);
	mpz_init(R[1]);
	mpz_init(R1[0]);
	mpz_init(R1[1]);

	mpz_init(k1);
	mpz_init(k2);
	mpz_init(b1);


	mpz_set(P[0], p1[0]);
	mpz_set(P[1], p1[1]);


	mpz_set(k2, *k);

	while (mpz_sgn(k2) == 1) {

		i = i + 1;


		b = mpz_tdiv_qr_ui(k1, b1, k2, 4);


		if (b == (long) 0) {
			b = mpz_tdiv_qr_ui(k1, b1, k2, 2);
			mpz_set(k2, k1);
			add_point(curv, P, P, Q, p);
			mpz_set(P[0],Q[0]);
			mpz_set(P[1],Q[1]);
		}


		else if (b == (long) 2) {
			b = mpz_tdiv_qr_ui(k1, b1, k2, 2);
			add_point(curv, P, P, Q, p);
			mpz_set(P[0],Q[0]);
			mpz_set(P[1],Q[1]);
			mpz_set(k2, k1);  			
		}	

		else if (b == (long) 1) {
			add_point(curv, Y, P, Q, p);
			mpz_set(Y[0],Q[0]);
			mpz_set(Y[1],Q[1]);					
			b = mpz_tdiv_qr_ui(k1, b1, k2, 2);
			mpz_set(k2, k1);
			add_point(curv, P, P, Q, p);
			mpz_set(P[0],Q[0]);
			mpz_set(P[1],Q[1]);
		
		}
				
		else if (b == (long) 3) {
		
			mpz_sub(P[1], *p, P[1]);
			add_point(curv, Y, P, Q, p);
			mpz_set(Y[0],Q[0]);
			mpz_set(Y[1],Q[1]);
			mpz_sub(P[1], *p, P[1]);
			mpz_add(k2,k2,one);
			b = mpz_tdiv_qr_ui(k1, b1, k2, 2);
			mpz_set(k2, k1);
			add_point(curv, P, P, Q, p);
			mpz_set(P[0],Q[0]);
			mpz_set(P[1],Q[1]);
		
		}	
		
	}

	mpz_set(result[0], Y[0]);
	mpz_set(result[1], Y[1]);


	mpz_clear(Y[0]);
	mpz_clear(Y[1]);
	mpz_clear(P[0]);
	mpz_clear(P[1]);
	mpz_clear(R[0]);
	mpz_clear(R[1]);
	mpz_clear(R1[0]);
	mpz_clear(R1[1]);

	mpz_clear(one);
	mpz_clear(k1);
	mpz_clear(k2);
	mpz_clear(b1);
}


/* implement scalar multiplication of the point p1 and k modulo p
   in the elliptic curve curv using the Window  method of size 4 256bit can deal*/
void point_mult_window(mpz_t * curv, mpz_t * p1, mpz_t * k1, mpz_t * result,
				mpz_t * p)
	
{
	mpz_t n,temp_t,ph[16][2],Y[2],p_temp[2];
	int *k,i,j,pow2;
	mpz_init(n);
	for(j=0 ; j<16 ;j++){
		mpz_init(ph[j][0]);
		mpz_init(ph[j][1]);
	}
	mpz_set(ph[0][0],p1[0]);
	mpz_set(ph[0][1],p1[1]);
	mpz_init(temp_t); 	
	mpz_init_set(Y[0],ph[0][0]);
	mpz_init_set(Y[1],ph[0][0]);
	mpz_init(p_temp[0]);
	mpz_init(p_temp[1]);
	mpz_set(n,*k1);
	
	i=-1;
	pow2=16;
	k=malloc(sizeof(int));
	
	while(mpz_sgn(n)>0)
	{
		i++;
	
		if(mpz_even_p(n)>0)
		{
			k[i]=0;
		}
	
		else
		{
			k[i]=mpz_fdiv_r_ui(temp_t,n,pow2);
			mpz_sub(n,n,temp_t);
		}

		mpz_fdiv_q_ui(temp_t,n,2);
		mpz_set(n,temp_t);
		k=realloc(k,(i+2)*sizeof(int));
	}
	for(j=0; j<15 ; j++){
			add_point(curv,ph[0],ph[j],ph[j+1],p);
	}
	if(k[i]){
		mpz_set(Y[0],ph[k[i]-1][0]);
		mpz_set(Y[1],ph[k[i]-1][1]);
	}

	for(j=i-1;j>=0;j--)
	{

		add_point(curv, Y, Y, p_temp, p);
		mpz_set(Y[0],p_temp[0]);
		mpz_set(Y[1],p_temp[1]);

		if(k[j]==1)
		{
			add_point(curv, Y, ph[0], p_temp, p);
			mpz_set(Y[0],p_temp[0]);
			mpz_set(Y[1],p_temp[1]);
		}
		
		else if(k[j]==3)
		{
			add_point(curv, Y, ph[2], p_temp, p);
			mpz_set(Y[0],p_temp[0]);
			mpz_set(Y[1],p_temp[1]);
		}
		
		else if(k[j]==5)
		{	
			add_point(curv, Y, ph[4], p_temp, p);
			mpz_set(Y[0],p_temp[0]);
			mpz_set(Y[1],p_temp[1]);
		}

		else if(k[j]==7)
		{
			add_point(curv, Y,ph[6], p_temp, p);
			mpz_set(Y[0],p_temp[0]);
			mpz_set(Y[1],p_temp[1]);
		}
		else if(k[j]==9)
		{
			add_point(curv, Y, ph[8], p_temp, p);
			mpz_set(Y[0],p_temp[0]);
			mpz_set(Y[1],p_temp[1]);
		}
		else if(k[j]==11)
		{
			add_point(curv, Y, ph[10], p_temp, p);
			mpz_set(Y[0],p_temp[0]);
			mpz_set(Y[1],p_temp[1]);
		}
		else if(k[j]==13)
		{
			add_point(curv, Y, ph[12], p_temp, p);
			mpz_set(Y[0],p_temp[0]);
			mpz_set(Y[1],p_temp[1]);
		}
		else if(k[j]==15)
		{
			add_point(curv, Y,ph[14], p_temp, p);
			mpz_set(Y[0],p_temp[0]);
			mpz_set(Y[1],p_temp[1]);
		}

	}	

	mpz_set(result[0],Y[0]);
	mpz_set(result[1],Y[1]);

		
	mpz_clear(n);
	mpz_clear(temp_t);
	mpz_clear(Y[0]);
	mpz_clear(Y[1]);	
	mpz_clear(p_temp[0]);
	mpz_clear(p_temp[1]);	
	for(j=0 ; j<16 ;j++){
		mpz_clear(ph[j][0]);
		mpz_clear(ph[j][1]);
	}
	free(k);
}

/* create a public key from the private key and the base point of the curve */
void my_public_key(mpz_t * curv, mpz_t * public1, mpz_t * basepoint,
				   mpz_t * private_key, mpz_t * p)
{

	point_mult(curv, basepoint, private_key, public1, p);

}


/* multiply the public key of one party with the private key of the other and
   generate the shared key */
void DH_shared_key(mpz_t * curv, mpz_t * shared_key, mpz_t * their_public,
				   mpz_t * private_key, mpz_t * p)
{
	point_mult(curv, their_public, private_key, shared_key, p);
}




/* create a public and a private key using the base point base */
void create_priv_and_public(mpz_t * curv, mpz_t * p, mpz_t * base,
					mpz_t * my_private, mpz_t * my_public)
{
	mpz_t k;

	mpz_init(k);

	myprimegenerator(&k);

	mpz_mod(*my_private, k, *p);

	my_public_key(curv, my_public, base, my_private, p);

	mpz_clear(k);
}


/* create a base point of order n and return the integers n, h, where nh = m (order of the curve) */
void domain_parameters(mpz_t * curv, mpz_t * base_point, mpz_t * p, mpz_t * m,
					   mpz_t * n, mpz_t * h)
{
	int f;
	mpz_t res, res1, cof;
	mpz_t p1[2];

	mpz_init(res);
	mpz_init(res1);
	mpz_init(cof);
	mpz_init(p1[0]);
	mpz_init(p1[1]);

	rand_point(curv, p, p1);

	mpz_set(res1, *m);

	f = mypollardrho(&res1, &res, &cof);

	while (f == 1) {
		mpz_set(res1, cof);
		f = mypollardrho(&res1, &res, &cof);
	}

	if (f == 0) {
		printf("could not find a factor, exit\n");
		exit(0);
	}

	mpz_tdiv_q(res, *m, res1);

	mpz_set(*n, res1);
	mpz_set(*h, res);

	point_mult(curv, p1, h, base_point, p);

	mpz_clear(res);
	mpz_clear(res1);
	mpz_clear(cof);
	mpz_clear(p1[0]);
	mpz_clear(p1[1]);

}

