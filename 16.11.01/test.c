#include "ECC.h"


void test_div_mul()
{
	mpz_t a, b, c, d, q, r;
	gmp_randstate_t state;

	mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d); mpz_init(q); mpz_init(r);
	gmp_randinit_default(state);


	//a ����
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
	//test_div_mul();
	mpz_t a, b, n;
	gmp_randstate_t state;

	mpz_init(a); mpz_init(b); mpz_init(n);
	gmp_randinit_default(state);

	mpz_urandomb(a, state, 448);
	mpz_urandomb(n, state, 224);
	mpz_set(b, a);

	START_WATCH;
	mpz_mod(a, a, n);
	STOP_WATCH;
	PRINT_TIME("mod ");
	gmp_printf("a = %Zx\n", a);

	START_WATCH;
	GFP_fast_reduction_p224(b, b, n);
	STOP_WATCH;
	PRINT_TIME("GFP ");
	gmp_printf("b = %Zx\n", b);

	mpz_clear(a); mpz_clear(n);
	gmp_randclear(state);

	return 0;
}