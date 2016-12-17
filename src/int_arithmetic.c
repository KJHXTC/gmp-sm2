// --------------------------------------------------------------------
//
//  File:        int_arithmetic.c
//  Date:        11/03
//  Last update: 11/03
//  Description: Integer arithmetic
//
//  (C) 2003, Elisavet Konstantinou & Yiannis Stamatiu & Christos Zaroliagis
//                 {konstane,stamatiu,zaro}@ceid.upatras.gr
//
// --------------------------------------------------------------------



#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "poly_arithmetic.h"
#include "int_arithmetic.h"


/* addition (x+y) modulo m */
void myzaddmod(mpz_t * res, mpz_t * x, mpz_t * y, mpz_t * m)
{
	mpz_t temp;

	mpz_init(temp);

	mpz_add(temp, *x, *y);
	mpz_mod(*res, temp, *m);

	mpz_clear(temp);
}


/* subtraction (x-y) modulo m */
void myzsubmod(mpz_t * res, mpz_t * x, mpz_t * y, mpz_t * m)
{
	mpz_t temp;

	mpz_init(temp);

	mpz_sub(temp, *x, *y);
	mpz_mod(*res, temp, *m);

	mpz_clear(temp);
}


/* multiplication (xy) modulo m */
void myzmulmod(mpz_t * res, mpz_t * x, mpz_t * y, mpz_t * m)
{
	mpz_t temp;

	mpz_init(temp);

	mpz_mul(temp, *x, *y);
	mpz_mod(*res, temp, *m);

	mpz_clear(temp);

}

void myzpowmod(mpz_t *res , mpz_t *x , unsigned long int exp ,mpz_t *m)
{
	mpz_t temp;
	mpz_init(temp);
	mpz_pow_ui(temp , *x ,exp);
	mpz_mod(*res , temp,*m);
	mpz_clear(temp);

}

/* multiplication (xy) (y unsigned long int) modulo m */
void myzsmulmod(mpz_t * res, mpz_t * x, unsigned long int y, mpz_t * m)
{
	mpz_t temp;

	mpz_init(temp);

	mpz_mul_ui(temp, *x, y);
	mpz_mod(*res, temp, *m);

	mpz_clear(temp);
}


/* division (x/y) modulo m */
void myzdivmod(mpz_t * res, mpz_t * x, mpz_t * y, mpz_t * m)
{
	int i;
	mpz_t h1;

	mpz_init(h1);

	i = mpz_invert(h1, *y, *m);

	if (i == 0) {
		printf("inverse is undefined!\n");
		exit(0);
	}

	else
		myzmulmod(res, x, &h1, m);

	mpz_clear(h1);
}


/* square root of a modulo p, r = sqrt(a) mod p */
void myzsqrtmod(mpz_t * r, mpz_t * a, mpz_t * p)
{
	long dP, rootl;
	mpz_t P1[3], helproot[2];

	mpz_init(P1[0]);
	mpz_init(P1[1]);
	mpz_init(P1[2]);
	mpz_init(helproot[0]);
	mpz_init(helproot[1]);

	mpz_set(P1[0], *a);
	mpz_neg(P1[0], P1[0]);
	mpz_set_ui(P1[1], 0);
	mpz_set_ui(P1[2], 1);

	dP = 2;

	rootl = 0;
	myRecurse(dP, p, P1, helproot, &rootl);

	mpz_set(*r, helproot[0]);

	mpz_clear(P1[0]);
	mpz_clear(P1[1]);
	mpz_clear(P1[2]);
	mpz_clear(helproot[0]);
	mpz_clear(helproot[1]);
}


/* generation of a prime number with size = (num_of_bits) bits
    (the definition of the global variable stat is necessary) */
void myprimegenerator(mpz_t * prime)
{
	mpz_t h1;

	mpz_init(h1);
	
	mpz_urandomb(h1, stat, bitlength);

	mpz_nextprime(*prime, h1);

	mpz_clear(h1);

}


/* function for the factorization of num to its two cofactors (res and cof)
   where res < cof. The function is based on the Pollard Rho algorithm.
   Returns 0 if it can't find a factorization, 1 if it finds one and
   2 if the number num is prime */
int mypollardrho1(mpz_t * num, mpz_t * res, mpz_t * cof)
{
	int primetest;

	mpz_t A, B;
	mpz_t h1, h2;

	mpz_init_set_ui(A, 1);
	mpz_init_set_ui(B, 2);
	mpz_init(h1);
	mpz_init(h2);

	primetest = mpz_probab_prime_p(*num, PRIMAL_TESTS);

	if (primetest != 0) {
		//   printf("number is prime\n");
		return 2;
	}

	else {
	  L1:
		primetest++;

		if (primetest > 1000)
			return 0;

		/*step 1 */
		mpz_mul(h1, A, A);
		mpz_mul(h1, h1, A);
		mpz_add_ui(h1, h1, 1);
		mpz_mod(A, h1, *num);

		/*step 2 */
		mpz_mul(h1, B, B);
		mpz_mul(h1, h1, B);
		mpz_add_ui(h1, h1, 1);
		mpz_mul(B, h1, h1);
		mpz_mul(B, B, h1);
		mpz_add_ui(h1, B, 1);
		mpz_mod(B, h1, *num);

		/*step 3 */
		mpz_sub(h1, B, A);
		mpz_gcd(h2, h1, *num);

		if (mpz_cmp_ui(h2, 1) == 0)
			goto L1;

		else if (mpz_cmp(h2, *num) == 0)
			return 0;

		else {
			mpz_div(h1, *num, h2);
			if (mpz_cmp(h1, h2) > 0) {
				mpz_set(*res, h2);
				mpz_set(*cof, h1);
			} else {
				mpz_set(*res, h1);
				mpz_set(*cof, h2);
			}
		}
	}

	mpz_clear(A);
	mpz_clear(B);
	mpz_clear(h1);
	mpz_clear(h2);

	return 1;
}


/* a slightly different algorithm from the above. In some cases, one of the two algorithms
   fails to find a factorization but the other can. */
int mypollardrho(mpz_t * num, mpz_t * res, mpz_t * cof)
{
	int f, primetest;

	mpz_t A, B;
	mpz_t h1, h2;

	mpz_init_set_ui(A, 1);
	mpz_init_set_ui(B, 2);
	mpz_init(h1);
	mpz_init(h2);

	primetest = mpz_probab_prime_p(*num, PRIMAL_TESTS);

	if (primetest != 0) {
		//  printf("number is prime\n");
		return 2;
	}

	else {
	  L1:
		primetest++;

		if (primetest > 2000)
			return 0;

		/*step 1 */
		mpz_mul(h1, A, A);
		mpz_add_ui(h1, h1, 1);
		mpz_mod(A, h1, *num);

		/*step 2 */
		mpz_mul(h1, B, B);
		mpz_add_ui(h1, h1, 1);
		mpz_mul(B, h1, h1);
		mpz_add_ui(h1, B, 1);
		mpz_mod(B, h1, *num);

		/*step 3 */
		mpz_sub(h1, B, A);
		mpz_gcd(h2, h1, *num);

		if (mpz_cmp_ui(h2, 1) == 0)
			goto L1;

		else if (mpz_cmp(h2, *num) == 0) {
			f = mypollardrho1(num, res, cof);
			return f;
		}

		else {
			mpz_div(h1, *num, h2);
			if (mpz_cmp(h1, h2) > 0) {
				mpz_set(*res, h2);
				mpz_set(*cof, h1);
			} else {
				mpz_set(*res, h1);
				mpz_set(*cof, h2);
			}
		}
	}

	mpz_clear(A);
	mpz_clear(B);
	mpz_clear(h1);
	mpz_clear(h2);

	return 1;
}

