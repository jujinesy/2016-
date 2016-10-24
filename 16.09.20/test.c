#include "gmp.h"
#include <stdlib.h>

int main()
{
	mpz_t a, b, c;
	
	mpz_init(a);
	//mpz_clear(a);
	a->_mp_d[0] = 0;
	a->_mp_d[1] = 0x00abcdef;

	mpz_init(b);
	//mpz_clear(b);
	b->_mp_d[0] = 0;
	b->_mp_d[1] = 0x12345678;

	mpz_init(c);


	mpz_mul(c, a, b);


	return 0;
}




// 1024 보다 작은 소수는 많다
// MR테스트 56번 한다