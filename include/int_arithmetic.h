//    This file is part of ecc-lib-2.0.
//
//    ecc-lib-2.0 is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    ecc-lib-2.0 is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ecc-lib-2.0; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// -----------------------------------------------------------------------------
//
//  File:        int_arithmetic.h
//  Date:        11/03
//  Last update: 11/03
//  Description: Integer arithmetic
//
//  (C) 2003, Elisavet Konstantinou & Yannis Stamatiou & Christos Zaroliagis
//                 {konstane,stamatiu,zaro}@ceid.upatras.gr
//
// -----------------------------------------------------------------------------


#ifndef INT_ARITHMETICH
#define INT_ARITHMETICH


#define PRIMAL_TESTS 10 /* number of tests for primality testing */
#define bitlength 256 /* size of the underlying finite field Fp */

extern gmp_randstate_t stat; /* global variable necessary for the generation of random integers */


/* addition (x+y) modulo m */
void myzaddmod(mpz_t *res, mpz_t *x, mpz_t *y, mpz_t *m);


/* subtraction (x-y) modulo m */
void myzsubmod(mpz_t *res, mpz_t *x, mpz_t *y, mpz_t *m);


/* multiplication (xy) modulo m */
void myzmulmod(mpz_t *res, mpz_t *x, mpz_t *y, mpz_t *m);


/* multiplication (xy) (y unsigned long int) modulo m */
void myzsmulmod(mpz_t *res, mpz_t *x, unsigned long int y, mpz_t *m);


/* division (x/y) modulo m */
void myzdivmod(mpz_t *res, mpz_t *x, mpz_t *y, mpz_t *m);


/* square root of a modulo p, r = sqrt(a) mod p */
void myzsqrtmod(mpz_t *r, mpz_t *a, mpz_t *p);


/* generation of a prime number with size = (num_of_bits) bits
    (the definition of the global variable stat is necessary) */
void myprimegenerator(mpz_t *prime);


/* function for the factorization of num to its two cofactors (res and cof)
   where res < cof. The function is based on the Pollard Rho algorithm.
   Returns 0 if it can't find a factorization, 1 if it finds one and
   2 if the number num is prime */
int mypollardrho1(mpz_t *num, mpz_t *res, mpz_t *cof);


/* a slightly different algorithm from the above. In some cases, one of the two algorithms
   fails to find a factorization but the other can. */
int mypollardrho(mpz_t *num, mpz_t *res, mpz_t *cof);


void myzpowmod(mpz_t *res , mpz_t *x , unsigned long int exp ,mpz_t *m);

#endif

