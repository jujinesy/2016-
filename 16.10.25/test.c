#include "ECC.h"


void test_div_mul()
{
	mpz_t a, b, c, d, q, r;
	gmp_randstate_t state;

	mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d); mpz_init(q); mpz_init(r);
	gmp_randinit_default(state);


	//a »ý¼º
	mpz_urandomb(a, state, 1024);
	mpz_urandomb(b, state, 1024);
	mpz_urandomb(d, state, 1024);

	START_WATCH;
	for (int i = 0; i < 100000; i++) {
		mpz_mul(c, a, b);
	}
	STOP_WATCH;
	PRINT_TIME("1024bit mul ");

	START_WATCH;
	for (int i = 0; i < 100000; i++) {
	//	mpz_urandomb(a, state, 1024);
	//	mpz_urandomb(b, state, 1024);
	//	mpz_urandomb(d, state, 1024); 
	//	
	//	mpz_mul(c, a, b);
		mpz_mod(r, c, d);
	}
	STOP_WATCH;
	PRINT_TIME("1024bit div ");
	//mpz_tdiv_qr(q, r, c, d);


	mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d); mpz_clear(q); mpz_clear(r);
	gmp_randclear(state);
}

int main()
{
	test_div_mul();
	return 0;
}