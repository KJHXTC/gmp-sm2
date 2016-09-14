// --------------------------------------------------------------------
//
//  File:        poly_arithmetic.c
//  Date:        11/03
//  Last update: 11/03
//  Description: Polynomial arithmetic
//
//  (C) 2003, Elisavet Konstantinou & Yiannis Stamatiu & Christos Zaroliagis
//                 {konstane,stamatiu,zaro}@ceid.upatras.gr
//
// --------------------------------------------------------------------



#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "int_arithmetic.h"
#include "poly_arithmetic.h"


void zpoly_copy(long da, mpz_t * za, mpz_t * zb, long *db)
{

	long i;



	*db = da;

	for (i = 0; i <= da; i++)
		mpz_set(zb[i], za[i]);

}

void zpoly_mul(long m, long n, mpz_t * za, mpz_t * zb, mpz_t * zc, long *p)
{

	long i, j, k;

	mpz_t zd, zai, zbk, zsum, zterm;

	mpz_init(zd);
	mpz_init(zai);
	mpz_init(zbk);
	mpz_init(zsum);
	mpz_init(zterm);

	*p = m + n;

	for (k = 0; k <= *p; k++) {

		mpz_set_ui(zsum, 0);

		for (i = 0; i <= k; i++) {

			j = k - i;

			if (i > m)
				mpz_set_ui(zai, 0);

			else
				mpz_set(zai, za[i]);

			if (j > n)
				mpz_set_ui(zbk, 0);

			else
				mpz_set(zbk, zb[j]);

			mpz_mul(zterm, zai, zbk);

			mpz_set(zd, zsum);

			mpz_add(zsum, zterm, zd);

		}

		mpz_set(zc[k], zsum);

	}

	mpz_clear(zd);
	mpz_clear(zai);
	mpz_clear(zbk);
	mpz_clear(zsum);
	mpz_clear(zterm);
}



void zpoly_div(long m, long n, mpz_t * zu, mpz_t * zv, mpz_t * zq, mpz_t * zr,
			   long *p, long *s)
{

	long j, jk, k, nk;

	mpz_t za, zb, zvn;

	mpz_init(za);
	mpz_init(zb);
	mpz_init(zvn);

	mpz_set(zvn, zv[n]);

	for (j = 0; j <= m; j++)

		mpz_set(zr[j], zu[j]);

	if (m < n) {

		*p = 0, *s = m;

		mpz_set_ui(zq[0], 0);

	}

	else {

		*p = m - n, *s = n - 1;

		for (k = *p; k >= 0; k--) {

			nk = n + k;

			mpz_pow_ui(za, zvn, k);

			mpz_mul(zq[k], zr[nk], za);

			for (j = nk - 1; j >= 0; j--) {

				jk = j - k;

				if (jk >= 0) {

					mpz_mul(za, zvn, zr[j]);

					mpz_mul(zb, zr[nk], zv[jk]);

					mpz_sub(zr[j], za, zb);

				}

				else {

					mpz_set(za, zr[j]);

					mpz_mul(zr[j], zvn, za);

				}

			}

		}

		while (*p > 0 && mpz_cmp_ui(zq[*p], 0l) == 0)
			*p = *p - 1;

		while (*s > 0 && mpz_cmp_ui(zr[*s], 0l) == 0)
			*s = *s - 1;

	}

	mpz_clear(za);
	mpz_clear(zb);
	mpz_clear(zvn);
}

void zpoly_mod(mpz_t * zp, mpz_t * za, long *da)
{

	long i;

	mpz_t zb;

	mpz_init(zb);

	for (i = 0; i <= *da; i++) {

		mpz_mod(zb, za[i], *zp);

		mpz_set(za[i], zb);

	}

	while (*da > 0 && mpz_cmp_ui(za[*da], 0l) == 0)
		*da = *da - 1;

	mpz_clear(zb);

}



void zpoly_pow(long degreeA, long degreem, mpz_t * zn, mpz_t * zp, mpz_t * zA,
			   mpz_t * zm, mpz_t * zs, long *ds)
{

	long dP, dq, dx = degreeA, i;

	mpz_t za, zb, zP[POLY_SIZE], zq[POLY_SIZE], zx[POLY_SIZE], zy[POLY_SIZE];

	mpz_init(za);
	mpz_init(zb);

	for (i = 0; i < POLY_SIZE; i++) {
		mpz_init(zP[i]);
		mpz_init(zq[i]);
		mpz_init(zx[i]);
		mpz_init(zy[i]);
	}

	*ds = 0;

	mpz_set(za, *zn);

	mpz_set_ui(zs[0], 1l);

	for (i = 0; i <= dx; i++)
		mpz_set(zx[i], zA[i]);

	while (mpz_cmp_ui(za, 0l) > 0) {

		if (mpz_tdiv_ui(za, 2l) == 1) {

			// s = (s * x) % m;

			zpoly_mul(*ds, dx, zs, zx, zP, &dP);

			zpoly_div(dP, degreem, zP, zm, zq, zs, &dq, ds);

			zpoly_mod(zp, zs, ds);

		}

		mpz_set(zb, za);

		mpz_tdiv_q_2exp(za, zb, 1l);

		if (mpz_cmp_ui(za, 0l) > 0) {

			// x = (x * x) % m;

			for (i = 0; i <= dx; i++)
				mpz_set(zy[i], zx[i]);

			zpoly_mul(dx, dx, zx, zy, zP, &dP);

			zpoly_div(dP, degreem, zP, zm, zq, zx, &dq, &dx);

			zpoly_mod(zp, zx, &dx);

		}

	}

	mpz_clear(za);
	mpz_clear(zb);

	for (i = 0; i < POLY_SIZE; i++) {
		mpz_clear(zP[i]);
		mpz_clear(zq[i]);
		mpz_clear(zx[i]);
		mpz_clear(zy[i]);
	}


}



void zpoly_sub(long da, long db, mpz_t * za, mpz_t * zb, mpz_t * zc, long *dc)
{

	long i;

	mpz_t zz;

	mpz_init(zz);

	mpz_set_ui(zz, 0);

	if (da >= db) {

		for (i = 0; i <= db; i++)

			mpz_sub(zc[i], za[i], zb[i]);

		for (i = db + 1; i <= da; i++)

			mpz_set(zc[i], za[i]);

		*dc = da;

	}

	else {

		for (i = 0; i <= da; i++)

			mpz_sub(zc[i], za[i], zb[i]);

		for (i = da + 1; i <= db; i++)

			mpz_sub(zc[i], zz, zb[i]);

		*dc = db;

	}

	mpz_clear(zz);

}

void zpoly_gcd(long degreeA, long degreeB, mpz_t * zp, mpz_t * zA, mpz_t * zB,
			   mpz_t * za, long *da)
{

	int nonzero = 0, zero;

	long db, dq, dr, i;

	mpz_t zc;

	mpz_t zb[POLY_SIZE], zq[POLY_SIZE], zr[POLY_SIZE];


	mpz_init(zc);

	for (i = 0; i < POLY_SIZE; i++) {
		mpz_init(zb[i]);
		mpz_init(zq[i]);
		mpz_init(zr[i]);
	}


	if (degreeA > degreeB) {

		*da = degreeA;

		db = degreeB;

		for (i = 0; i <= *da; i++)
			mpz_set(za[i], zA[i]);

		for (i = 0; i <= db; i++)
			mpz_set(zb[i], zB[i]);

	}

	else {

		*da = degreeB;

		db = degreeA;

		for (i = 0; i <= *da; i++)
			mpz_set(za[i], zB[i]);

		for (i = 0; i <= db; i++)
			mpz_set(zb[i], zA[i]);

	}

	for (i = 0; i <= db && !nonzero; i++)

		nonzero = mpz_cmp_ui(zb[i], 0l) != 0;

	while (nonzero) {

		zpoly_div(*da, db, za, zb, zq, zr, &dq, &dr);

		for (i = 0; i <= dr; i++) {

			mpz_set(zc, zr[i]);

			mpz_mod(zr[i], zc, *zp);

		}

		zero = 1;

		for (i = dr; i >= 0 && zero; i--) {

			zero = mpz_cmp_ui(zr[i], 0l) == 0;

			if (zero && dr > 0)
				dr--;

		}

		for (i = 0; i <= db; i++)
			mpz_set(za[i], zb[i]);

		*da = db;

		for (i = 0; i <= dr; i++)
			mpz_set(zb[i], zr[i]);

		db = dr;

		nonzero = 0;

		for (i = 0; i <= db && !nonzero; i++)

			nonzero = mpz_cmp_ui(zb[i], 0l) != 0;

	}


	mpz_clear(zc);

	for (i = 0; i < POLY_SIZE; i++) {
		mpz_clear(zb[i]);
		mpz_clear(zq[i]);
		mpz_clear(zr[i]);
	}

}



void zpoly_print(long da, mpz_t * za)
{

	long i;

	for (i = da; i >= 0; i--) {

		mpz_out_str(stdout, 10, za[i]);

		printf(" ");

	}

	printf("\n");

}



void zpoly_ext_euclid(long dg, long dh, mpz_t * zp, mpz_t * zg, mpz_t * zh,
					  mpz_t * zs, mpz_t * zt, mpz_t * zd, long *ds, long *dt,
					  long *dd)
{

	long da, dq, dr, ds1 = 0, ds2 = 0, dt1 = 0, dt2 = 0, i;

	mpz_t za[POLY_SIZE], zb[POLY_SIZE];

	mpz_t zq[POLY_SIZE], zr[POLY_SIZE];

	mpz_t zs1[POLY_SIZE], zs2[POLY_SIZE];

	mpz_t zt1[POLY_SIZE], zt2[POLY_SIZE];



	for (i = 0; i < POLY_SIZE; i++) {

		mpz_init(za[i]);
		mpz_init(zb[i]);
		mpz_init(zq[i]);
		mpz_init(zr[i]);

		mpz_init(zs1[i]);
		mpz_init(zs2[i]);
		mpz_init(zt1[i]);
		mpz_init(zt2[i]);

	}


	if (dh == 0 && mpz_cmp_ui(zh[0], 0l) == 0) {

		zpoly_copy(dg, zg, zd, dd);

		*ds = *dt = 0;

		mpz_set_ui(zs[0], 1l);

		mpz_set_ui(zt[0], 0l);

	}



	mpz_set_ui(zs2[0], 1l);

	mpz_set_ui(zs1[0], 0l);

	mpz_set_ui(zt2[0], 0l);

	mpz_set_ui(zt1[0], 1l);

	while (dh != 0 || mpz_cmp_ui(zh[0], 0l) != 0) {

		zpoly_div(dg, dh, zg, zh, zq, zr, &dq, &dr);

		zpoly_mod(zp, zq, &dq);

		zpoly_mod(zp, zr, &dr);

		zpoly_mul(dq, ds1, zq, zs1, za, &da);

		zpoly_sub(ds2, da, zs2, za, zs, ds);

		zpoly_mul(dq, dt1, zq, zt1, za, &da);

		zpoly_sub(dt2, da, zt2, za, zt, dt);

		zpoly_mod(zp, zs, ds);

		zpoly_mod(zp, zt, dt);

		zpoly_copy(dh, zh, zg, &dg);

		zpoly_copy(dr, zr, zh, &dh);

		zpoly_copy(ds1, zs1, zs2, &ds2);

		zpoly_copy(*ds, zs, zs1, &ds1);

		zpoly_copy(dt1, zt1, zt2, &dt2);

		zpoly_copy(*dt, zt, zt1, &dt1);

#ifdef DEBUG

		printf("q  = ");
		zpoly_print(dq, zq);

		printf("r  = ");
		zpoly_print(dr, zr);

		printf("s  = ");
		zpoly_print(*ds, zs);

		printf("t  = ");
		zpoly_print(*dt, zt);

		printf("g  = ");
		zpoly_print(dg, zg);

		printf("h  = ");
		zpoly_print(dh, zh);

		printf("s2 = ");
		zpoly_print(ds2, zs2);

		printf("s1 = ");
		zpoly_print(ds1, zs1);

		printf("t2 = ");
		zpoly_print(dt2, zt2);

		printf("t1 = ");
		zpoly_print(dt1, zt1);

#endif

	}

	zpoly_copy(dg, zg, zd, dd);

	zpoly_copy(ds2, zs2, zs, ds);

	zpoly_copy(dt2, zt2, zt, dt);


	for (i = 0; i < POLY_SIZE; i++) {

		mpz_clear(za[i]);
		mpz_clear(zb[i]);
		mpz_clear(zq[i]);
		mpz_clear(zr[i]);

		mpz_clear(zs1[i]);
		mpz_clear(zs2[i]);
		mpz_clear(zt1[i]);
		mpz_clear(zt2[i]);

	}

}


/* find the roots of the polynomial zA modulo zp */
void Recurse(long degreeA, mpz_t * zp, mpz_t * zA, mpz_t * zroot,
			 long *rootSize)
{

	long dd, degreeB, dq, dr, du = 1, i, flag = 0;

	mpz_t zD, za, zb, zc, ze;

	mpz_t zn, x0;

	mpz_t zB[POLY_SIZE], zd[POLY_SIZE];

	mpz_t zq[POLY_SIZE], zr[POLY_SIZE];

	mpz_t zu[2];


	mpz_init(zD);
	mpz_init(za);
	mpz_init(zb);
	mpz_init(zc);
	mpz_init(ze);
	mpz_init(zn);
	mpz_init(x0);

	for (i = 0; i < POLY_SIZE; i++) {
		mpz_init(zB[i]);
		mpz_init(zd[i]);
		mpz_init(zq[i]);
		mpz_init(zr[i]);
	}

	mpz_init(zu[0]);
	mpz_init(zu[1]);


	mpz_set_ui(x0, (long) 102);


	if (mpz_cmp_ui(zA[degreeA], (long) 1) != 0) {
		for (i = 0; i < (degreeA + 1); i++) {
			myzdivmod(&zA[i], &zA[i], &zA[degreeA], zp);

		}
	}


	while (degreeA != 1) {
		i = 0;
		do {

			mpz_sub_ui(za, *zp, 1l);
			mpz_tdiv_q_2exp(zn, za, 1l);

			mpz_add_ui(x0, x0, (long) 1);

			mpz_set(zu[0], x0);
			mpz_set_ui(zu[1], 1l);

			zpoly_mod(zp, zu, &du);

			zpoly_pow(du, degreeA, &zn, zp, zu, zA, zd, &dd);
			zpoly_mod(zp, zd, &dd);

			mpz_sub_ui(zd[0], zd[0], 1l);

			zpoly_gcd(dd, degreeA, zp, zd, zA, zB, &degreeB);

			zpoly_mod(zp, zB, &degreeB);


			if ((mpz_cmp_ui(x0, (long) 200) == 1)
				|| *rootSize > (degreeA - 1)) {
				flag = 1;

				goto L1;
			}


		} while ((degreeB == 0 || degreeB == degreeA));


		if (degreeB >= 1 && flag != 1) {

			Recurse(degreeB, zp, zB, zroot, rootSize);

			zpoly_div(degreeA, degreeB, zA, zB, zq, zr, &dq, &dr);
			zpoly_mod(zp, zq, &dq);
			zpoly_copy(dq, zq, zA, &degreeA);

			Recurse(degreeA, zp, zA, zroot, rootSize);

		}

	}


	if (degreeA == 1) {
		mpz_invert(za, zA[1], *zp);
		mpz_mul(zb, zA[0], za);
		mpz_neg(zb, zb);
		mpz_mod(zroot[*rootSize], zb, *zp);

		*rootSize = *rootSize + 1;

	}



  L1:

	mpz_clear(zD);
	mpz_clear(za);
	mpz_clear(zb);
	mpz_clear(zc);
	mpz_clear(ze);
	mpz_clear(zn);
	mpz_clear(x0);

	for (i = 0; i < POLY_SIZE; i++) {
		mpz_clear(zB[i]);
		mpz_clear(zd[i]);
		mpz_clear(zq[i]);
		mpz_clear(zr[i]);
	}

	mpz_clear(zu[0]);
	mpz_clear(zu[1]);

}


/* find a single root of polynomial zA modulo zp */
void myRecurse(long degreeA, mpz_t * zp, mpz_t * zA, mpz_t * zroot,
			   long *rootSize)
{

	long dd, degreeB, dq, dr, du = 1, i, flag = 0;

	mpz_t zn, x0, za;

	mpz_t zB[POLY_SIZE], zd[POLY_SIZE];

	mpz_t zq[POLY_SIZE], zr[POLY_SIZE];

	mpz_t zu[2];


	mpz_init(zn);
	mpz_init(x0);
	mpz_init(za);

	for (i = 0; i < POLY_SIZE; i++) {
		mpz_init(zB[i]);
		mpz_init(zd[i]);
		mpz_init(zq[i]);
		mpz_init(zr[i]);
	}

	mpz_init(zu[0]);
	mpz_init(zu[1]);


	mpz_set_ui(x0, (long) 102);


	if (mpz_cmp_ui(zA[degreeA], (long) 1) != 0) {
		for (i = 0; i < (degreeA + 1); i++) {
			myzdivmod(&zA[i], &zA[i], &zA[degreeA], zp);

		}
	}

	while (degreeA != 1) {
		i = 0;
		do {

			mpz_sub_ui(za, *zp, 1l);
			mpz_tdiv_q_2exp(zn, za, 1l);

			mpz_add_ui(x0, x0, (long) 1);

			mpz_set(zu[0], x0);
			mpz_set_ui(zu[1], 1l);

			i++;

			zpoly_mod(zp, zu, &du);

			zpoly_pow(du, degreeA, &zn, zp, zu, zA, zd, &dd);
			zpoly_mod(zp, zd, &dd);

			mpz_sub_ui(zd[0], zd[0], 1l);

			zpoly_gcd(dd, degreeA, zp, zd, zA, zB, &degreeB);

			zpoly_mod(zp, zB, &degreeB);


			if (i > 50 || *rootSize == 1) {
				flag = 1;

				goto L1;
			}


		} while ((degreeB == 0 || degreeB == degreeA));


		if (degreeB >= 1 && flag != 1) {
			if (degreeB > degreeA / 2) {
				zpoly_div(degreeA, degreeB, zA, zB, zq, zr, &dq, &dr);
				zpoly_mod(zp, zq, &dq);
				zpoly_copy(dq, zq, zA, &degreeA);

				Recurse(degreeA, zp, zA, zroot, rootSize);
			}

			else {
				zpoly_copy(degreeB, zB, zA, &degreeA);

				Recurse(degreeA, zp, zA, zroot, rootSize);

			}
		}

	}


	if (degreeA == 1) {
		mpz_invert(za, zA[1], *zp);
		mpz_mul(x0, zA[0], za);
		mpz_neg(x0, x0);
		mpz_mod(zroot[*rootSize], x0, *zp);

		*rootSize = *rootSize + 1;

	}

  L1:

	mpz_clear(zn);
	mpz_clear(x0);
	mpz_clear(za);

	for (i = 0; i < POLY_SIZE; i++) {
		mpz_clear(zB[i]);
		mpz_clear(zd[i]);
		mpz_clear(zq[i]);
		mpz_clear(zr[i]);
	}

	mpz_clear(zu[0]);
	mpz_clear(zu[1]);
}

