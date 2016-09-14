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
//  File:        poly_arithmetic.h
//  Date:        11/03
//  Last update: 11/03
//  Description: Polynomial arithmetic
//
//  (C) 2003, Elisavet Konstantinou & Yannis Stamatiou & Christos Zaroliagis
//                 {konstane,stamatiu,zaro}@ceid.upatras.gr
//
// -----------------------------------------------------------------------------


#ifndef POLY_ARITHMETICH
#define POLY_ARITHMETICH

#define POLY_SIZE 181


void zpoly_copy(long da, mpz_t *za, mpz_t *zb, long *db);


void zpoly_mul(long m, long n, mpz_t *za, mpz_t *zb, mpz_t *zc, long *p);


void zpoly_div(long m, long n, mpz_t *zu, mpz_t *zv, mpz_t *zq, mpz_t *zr, long *p, long *s);


void zpoly_mod(mpz_t *zp, mpz_t *za, long *da);


void zpoly_pow(long degreeA, long degreem, mpz_t *zn, mpz_t *zp, mpz_t *zA, mpz_t *zm, mpz_t *zs, long *ds);


void zpoly_sub(long da, long db, mpz_t *za, mpz_t *zb, mpz_t *zc, long *dc);


void zpoly_gcd(long degreeA, long degreeB, mpz_t *zp, mpz_t *zA, mpz_t *zB, mpz_t *za, long *da);


void zpoly_print(long da, mpz_t *za);


void zpoly_ext_euclid(long dg, long dh, mpz_t *zp, mpz_t *zg, mpz_t *zh, mpz_t *zs,
                      mpz_t *zt, mpz_t *zd, long *ds, long *dt, long *dd);


/* find the roots of the polynomial zA modulo zp */
void Recurse(long degreeA, mpz_t *zp, mpz_t *zA, mpz_t *zroot, long *rootSize);


/* find a single root of polynomial zA modulo zp */
void myRecurse(long degreeA, mpz_t *zp, mpz_t *zA, mpz_t *zroot, long *rootSize);

  
#endif

