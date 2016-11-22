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
	mpz_urandomb(d, state, 1024);;

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

void GFP_fast_reduction_p224_test()
{
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
	for (int i = 0; i<10000; i++)
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
}

void ecc_affine_add_doubling_test()
{
	GFP_POINT p, q, r, pq, p2, r2;
	mpz_t prime, co_a;
	mpz_init(prime); mpz_init(co_a);



	GFP_point_init(&p);
	GFP_point_init(&q);
	GFP_point_init(&r);
	GFP_point_init(&r2);
	GFP_point_init(&pq);
	GFP_point_init(&p2);

	//p
	mpz_set_str(p.x, "b70e0cbd6bb4bf7f321390b94a03c1d356c21122343280d6115c1d21", 16);
	mpz_set_str(p.y, "bd376388b5f723fb4c22dfe6cd4375a05a07476444d5819985007e34", 16);
	p.point_at_infinity = 0;
	mpz_set_str(co_a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFE",16);
	mpz_set_str(prime, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001", 16);


	//q
	mpz_set_str(q.x, "3e54b1067d79dab2bd159a61f03e7bcb62851c5d3418a8bc5cad16aa", 16);
	mpz_set_str(q.y, "5986453b559d4ce6b24a4f10a87efe19eaf49a6db877bda6b791adb3", 16);
	q.point_at_infinity = 0;

	//addition
	GFP_affine_addition(&r, &p, &q, co_a, prime);

	gmp_printf("add rx = %Zx\n", r.x);
	gmp_printf("add ry = %Zx\n", r.y);


	//doubling
	GFP_affine_doubling(&r2, &p, co_a, prime);
	gmp_printf("dbl rx = %Zx\n", r2.x);
	gmp_printf("dbl ry = %Zx\n", r2.y);


	
	GFP_point_clear(&p);
	GFP_point_clear(&q);
	GFP_point_clear(&r);
	GFP_point_clear(&r2);
	GFP_point_clear(&pq);
	GFP_point_clear(&p2);
	mpz_clear(prime);
	mpz_clear(co_a);
}

int main()
{
	//test_div_mul();
	ecc_affine_add_doubling_test();

	return 0;
}


//be5f