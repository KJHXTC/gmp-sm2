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
//  File:        ec_operations.h
//  Date:        11/03
//  Last update: 04/10
//  Description: Basic operations on elliptic curves' group
//
//  (C) 2003, Elisavet Konstantinou & Yannis Stamatiou & Christos Zaroliagis
//                 {konstane,stamatiu,zaro}@ceid.upatras.gr
//
//   Extended by: Bhanu Prakash & Pratik Poddar & Bernard Menezes 
//                 {prakashb,pratik, bernard}@cse.iitb.ac.in
//
// -----------------------------------------------------------------------------


#ifndef EC_OPERATIONSH
#define EC_OPERATIONSH


# define MAX_CHAR 100 /* maximum number of characters in a string */
    
extern gmp_randstate_t stat;
 

/* generates randomly the a and b coefficients (curv[0] and curv[1]
   respectively) of the elliptic curve modulo the prime p */
void rand_curve (mpz_t *curv, mpz_t *p);
  

/* print the coefficients of an elliptic curve */ 
void print_curve( mpz_t *curv);


/* print a point */
void print_point( mpz_t *point);
 
 
/* addition of two points p1 and p2 modulo p that belong in the 
  elliptic curve curv. The result of the addition is the point p3 */
void add_point( mpz_t *curv, mpz_t *p1, mpz_t *p2, mpz_t *p3, mpz_t *p);


/* implement scalar multiplication of the point p1 and k modulo p 
   in the elliptic curve curv using the binary method */ 
void point_mult(mpz_t *curv, mpz_t *p1, mpz_t *k, mpz_t *result, mpz_t *p);


/* implement scalar multiplication of the point p1 and k modulo p
   in the elliptic curve curv using the naf method */
void point_mult_naf(mpz_t *curv, mpz_t *p1, mpz_t *k, mpz_t *result, mpz_t *p);


/* implement scalar multiplication of the point p1 and k modulo p
   in the elliptic curve curv using the Window  method of size 3 */
void point_mult_window(mpz_t *curv, mpz_t *p1, mpz_t *k1, mpz_t *result, mpz_t *p);


/* create a public key from the private key and the base point of the curve */
void my_public_key(mpz_t *curv, mpz_t *public1, mpz_t *basepoint, mpz_t *private_key, mpz_t *p);


/* multiply the public key of one party with the private key of the other and 
   generate the shared key */
void DH_shared_key(mpz_t *curv, mpz_t *shared_key, mpz_t *their_public, mpz_t *private_key, mpz_t *p);
 

/* generate a random point base in the elliptic curve curv modulo p */
void rand_point(mpz_t *curv, mpz_t *p, mpz_t *base);

 
/* create a public and a private key using the base point base */
void create_priv_and_public(mpz_t *curv, mpz_t *p, mpz_t *base, mpz_t *my_private, mpz_t *my_public);


/* create a base point of order n and return the integers n, h, where nh = m (order of the curve) */
void domain_parameters(mpz_t *curv, mpz_t *base_point, mpz_t *p, mpz_t *m, mpz_t *n, mpz_t *h);

#endif


