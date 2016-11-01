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
	//test_div_mul();
	mpz_t a, b, c, n;
	gmp_randstate_t state;

	mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(n);
	gmp_randinit_default(state);

	mpz_urandomb(a, state, 440);
	mpz_urandomb(n, state, 224);
	gmp_printf("a = %Zx\n", a);

	n->_mp_d[6] = 0xffffffff;
	n->_mp_d[5] = 0xffffffff;
	n->_mp_d[4] = 0xffffffff;
	n->_mp_d[3] = 0xffffffff;
	n->_mp_d[2] = 0x0;
	n->_mp_d[1] = 0x0;
	n->_mp_d[0] = 0x1;
	n->_mp_size = 7;

	START_WATCH;
	for(int i=0;i<10000;i++)
		mpz_mod(b, a, n);
	STOP_WATCH;
	PRINT_TIME("mod ");
	gmp_printf("b = %Zx\n", b);




	START_WATCH;
	for (int i = 0; i<10000; i++)
		GFP_fast_reduction_p224(c, a, n);
	STOP_WATCH;
	PRINT_TIME("GFP ");
	gmp_printf("c = %Zx\n", c);

	mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(n);
	gmp_randclear(state);

	return 0;
}