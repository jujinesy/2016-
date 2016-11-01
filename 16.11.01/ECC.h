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


// a mod n =r
int GFP_fast_reduction_p224(mpz_t r, const mpz_t a, const mpz_t n);
int GFP_fast_reduction_p256(mpz_t r, const mpz_t a, const mpz_t n);



// r = a+b mod n
void mpz_add_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n);
// r = a-b mod n
void mpz_sub_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n);
// r = a*b mod n
void mpz_mul_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n);