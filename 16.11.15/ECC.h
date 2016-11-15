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

typedef struct _NAF_RECORDING_
{
	char naf_scalar[1024];
	int	naf_len;
}NAF_RECORDING;





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
//R = P+Q
int GFP_affine_addition(GFP_POINT *r, GFP_POINT *p, GFP_POINT *q, const mpz_t coefficient_a, const mpz_t prime);



// R=kP LtoR binary
//int GFP_LtoR_binary(GFP_POINT *r, )
// R = kP LtoR NAF

// R=kP LtoR binary
int GFP_LtoR_binary(GFP_POINT *r, mpz_t k, GFP_POINT *p, const mpz_t coefficient_a, const mpz_t prime);

//R=kP LtoR NAF
int GFP_naf_recording(NAF_RECORDING *rk, mpz_t k);
int GFP_affine_addition(GFP_POINT *r, NAF_RECORDING *k, GFP_POINT *p, const mpz_t coefficient_a);