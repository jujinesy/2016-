#include "../include/gmp.h"
#include <stdlib.h>

int main()
{
	mpz_t p, q;
	gmp_randstate_t state;

	mpz_init(p);
	gmp_randinit_default(state);
	mpz_urandomb(p, state, 1024);
	p->_mp_d[0] |= 1;
	while (1)
	{
		if (mpz_probab_prime_p(p, 56))
			break;
		mpz_add_ui(p, p, 2);
	}
	gmp_printf("%Zd\n", p);
	mpz_clear(p);

	getchar();
	return 0;
}