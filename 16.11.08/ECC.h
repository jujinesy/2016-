//환경설정 PATH=$(SolutionDir)/include;%PATH%
#pragma comment(lib, "../include/gmp.lib") 
#include "../include/gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

clock_t elapsed;
float sec;

#define START_WATCH \
{\
 elapsed = -clock(); \
}\

#define STOP_WATCH \
{\
 elapsed += clock();\
 sec = (float)elapsed/CLOCKS_PER_SEC;\
}\

#define PRINT_TIME(qstr) \
{\
 printf("\n[%s: %.5f s]\n",qstr,sec);\
}


typedef struct _GFP_POINT_
{
	mpz_t x;
	mpz_t y;
	int point_at_infinity;
}GFP_POINT;



// a mod n =r
int GFP_fast_reduction_p224(mpz_t c, const mpz_t a, const mpz_t n);
int GFP_fast_reduction_p256(mpz_t c, const mpz_t a, const mpz_t n);



// r = a+b mod n
void mpz_add_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n);
// r = a-b mod n
void mpz_sub_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n);
// r = a*b mod n
void mpz_mul_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n);



//ECC computation
int GFP_point_init(GFP_POINT *p);
int GFP_point_clear(GFP_POINT *p);

//R = 2P
int GFP_affine_doubling(GFP_POINT *r, GFP_POINT *p, const mpz_t coefficient_a, const mpz_t prime);
//R = 
int GFP_affine_addition(GFP_POINT *r, GFP_POINT *p, GFP_POINT *q, const mpz_t coefficient_a, const mpz_t prime);


//R = 2P
void GFP_Affine_Doubling2(GFP_POINT *r, GFP_POINT *p, mpz_t a, const mpz_t prime);
//R = 
void GFP_Affine_Addition2(GFP_POINT *r, GFP_POINT *p, GFP_POINT *q, mpz_t a, const mpz_t prime);