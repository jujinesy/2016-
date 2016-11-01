#include "ECC.h"



// a mod n =r
int GFP_fast_reduction_p224(mpz_t c, const mpz_t a, const mpz_t n)
{
	//mpz_t s1;

	//mpz_init(s1);

	//s1->_mp_d[0] = a->_mp_d[0];
	//s1->_mp_d[1] = a->_mp_d[1];
	//s1->_mp_d[2] = a->_mp_d[2];
	//s1->_mp_size = 3;

	mpz_t s1, s2, s3, s4, s5;

	if (mpz_size(a) > 14) {
		mpz_mod(c, a, n);
		return;
	}

	mpz_init2(s1, 500);
	mpz_init2(s2, 256);
	mpz_init2(s3, 256);
	mpz_init2(s4, 256);
	mpz_init2(s5, 256);

	memset(s1->_mp_d, 0, 56);
	mpz_set(s1, a);
	s1->_mp_size = 7;

	s2->_mp_d[0] = s2->_mp_d[1] = s2->_mp_d[2] = 0;
	s3->_mp_d[0] = s3->_mp_d[1] = s3->_mp_d[2] = 0;

	s2->_mp_d[3] = s4->_mp_d[0] = s1->_mp_d[7];
	s2->_mp_d[4] = s4->_mp_d[1] = s1->_mp_d[8];
	s2->_mp_d[5] = s4->_mp_d[2] = s1->_mp_d[9];
	s2->_mp_d[6] = s4->_mp_d[3] = s1->_mp_d[10];
	s2->_mp_size = 7;

	s3->_mp_d[3] = s4->_mp_d[4] = s5->_mp_d[0] = s1->_mp_d[11];
	s3->_mp_d[4] = s4->_mp_d[5] = s5->_mp_d[1] = s1->_mp_d[12];
	s3->_mp_d[5] = s4->_mp_d[6] = s5->_mp_d[2] = s1->_mp_d[13];
	s3->_mp_size = 6;
	s4->_mp_size = 7;
	s5->_mp_size = 3;

	mpz_add(s2, s2, s1);
	mpz_add(s2, s2, s3);
	mpz_sub(s2, s2, s4);
	mpz_sub(s2, s2, s5);

	if (mpz_cmp(s2, n) >= 0) {
		while (mpz_cmp(s2, n) >= 0) {
			mpz_sub(s2, s2, n);
		}
	}
	else if (mpz_sgn(s2)<0) {
		while (mpz_sgn(s2)<0) {
			mpz_add(s2, s2, n);
		}
	}

	if (mpz_sgn(a)<0) mpz_sub(c, n, s2);
	else             mpz_set(c, s2);

	mpz_clear(s1); mpz_clear(s2); mpz_clear(s3); mpz_clear(s4); mpz_clear(s5);
	return 0;
}

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