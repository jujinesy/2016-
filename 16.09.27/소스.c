#include "gmp.h"
#include <stdlib.h>

typedef struct _RSA_PRIKEY_ {
	mpz_t n;
	mpz_t e;
}RSA_PUBKEY;

typedef struct _RSA_PRIKEY_2 {
	mpz_t n;
	mpz_t e;
	mpz_t p;
	mpz_t q;
	mpz_t d;
};

//void RSA_KEYGEN(RSA_PUBKEY * pubkey, RSA_PUBKEY *prikey )
int main()
{
	mpz_t p, q, n, e, d, tmp;
	gmp_randstate_t state;

	mpz_init(p);
	mpz_init(q);
	mpz_init(n);
	mpz_init(e);
	mpz_init(d);
	mpz_init(tmp);
	gmp_randinit_default(state);

	//e 대입  e=0x10001
	e->_mp_d[0] = 0x10001;
	e->_mp_size = 1;
	gmp_printf("e = %Zx\n", e);

	//mpz_set_ui(e, 10001); 와 동일

	mpz_urandomb(p, state, 1024);
	p->_mp_d[0] |= 1;
	// mpz_setbit(p, 0);

	p->_mp_d[31] |= 0x80000000;
	// mpz_setbit(p, 1023); // 쓰레기값 0으로 초기화도 해줌
	while (1)
	{
		if (mpz_probab_prime_p(p, 56)) {  //1확률적 소수 2는 100퍼 소수임 0은 소수아님 
			mpz_sub_ui(d, p, 1);
			mpz_gcd(n, e, d);

			//if (!mpz_cmp_si(n, 1))
			if ((n->_mp_size == 1) && (n->_mp_d[0] == 1))
				break;

			break;
		}
		mpz_add_ui(p, p, 2);
	}
	gmp_printf("p = %Zx\n", p);










	
	mpz_urandomb(q, state, 1024);
	q->_mp_d[0] |= 1;
	// mpz_setbit(p, 0);

	q->_mp_d[31] |= 0x80000000;
	// mpz_setbit(p, 1023); // 쓰레기값 0으로 초기화도 해줌
	while (1)
	{
		if (mpz_probab_prime_p(q, 56)) {  //1확률적 소수 2는 100퍼 소수임 0은 소수아님 
			mpz_sub_ui(d, q, 1);
			mpz_gcd(n, e, d);

			//if (!mpz_cmp_si(n, 1))
			if ((n->_mp_size == 1) && (n->_mp_d[0] == 1))
				break;

			break;
		}
		mpz_add_ui(q, q, 2);
	}
	gmp_printf("q = %Zx\n", q);


	mpz_mul(n, p, q);
	gmp_printf("n = %Zx\n", n);

	mpz_sub_ui(d, p, 1);
	mpz_sub_ui(tmp, q, 1);


	mpz_mul(tmp, tmp, d);
	mpz_invert(d, e, tmp);
	gmp_printf("d = %Zx\n", d);
	
	mpz_mul(e, e, d);
	mpz_mod(tmp, e, tmp);
	gmp_printf("e*d mod phi(n) = %Zx\n", tmp);

	mpz_clear(p); mpz_clear(q); mpz_clear(n); mpz_clear(e); mpz_clear(d); mpz_clear(tmp);
	gmp_randclear(state);

	getchar();
	return 0;
}