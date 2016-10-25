#include "ECC.h"


// r = a+b mod n
void mpz_add_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n)
{
	mpz_add(r, a, b);
	if (mpz_cmp(r, n) >= 0) mpz_sub(r, r, n);
	//mpz_mod(r, r, n);
}
// r = a-b mod n
void mpz_sub_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n)
{
	mpz_sub(r, a, b);
	if (mpz_sgn(r) < 0)	mpz_add(r, r, n);
	//mpz_mod(r, r, n);
}
// r = a*b mod n
void mpz_mul_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n)
{
	mpz_mul(r, a, b);
	mpz_mod(r, r, n);
}