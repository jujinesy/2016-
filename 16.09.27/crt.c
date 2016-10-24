#include "../include/gmp.h"
#include <stdlib.h>

//1. �ӵ�����
//2. Ű���� �߰�
//dp = d mod p-1
//dq = d mod q-1
//		-1
//qinv = q mod p
//
//3 ��ȣȭ
// e
//m mod n = c
//
//4 ��ȣȭ
// d
//c mod n = m
//
//
//crt
//            dp
//x = (c mod p) mod p
//
//             dq
//y = (c mod p) mod q
//
//m = y + (x-y)qinv*q mod n
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
	mpz_t p, q, n, e, d, tmp, dp, dq, qinv, m, c, x, y;
	gmp_randstate_t state;

	mpz_init(p); mpz_init(q); mpz_init(n); mpz_init(e); mpz_init(d); mpz_init(tmp); mpz_init(dp); mpz_init(dq); mpz_init(qinv); mpz_init(m); mpz_init(c); mpz_init(x); mpz_init(y);
	gmp_randinit_default(state);

	//e ����  e=0x10001
	e->_mp_d[0] = 0x10001;
	e->_mp_size = 1;
	gmp_printf("e = %Zx\n", e);

	//mpz_set_ui(e, 10001); �� ����

	mpz_urandomb(p, state, 1024);
	p->_mp_d[0] |= 1;
	// mpz_setbit(p, 0);

	p->_mp_d[31] |= 0x80000000;
	// mpz_setbit(p, 1023); // �����Ⱚ 0���� �ʱ�ȭ�� ����
	while (1)
	{
		if (mpz_probab_prime_p(p, 56)) {  //1Ȯ���� �Ҽ� 2�� 100�� �Ҽ��� 0�� �Ҽ��ƴ� 
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
	// mpz_setbit(p, 1023); // �����Ⱚ 0���� �ʱ�ȭ�� ����
	while (1)
	{
		if (mpz_probab_prime_p(q, 56)) {  //1Ȯ���� �Ҽ� 2�� 100�� �Ҽ��� 0�� �Ҽ��ƴ� 
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

	//mpz_mul(e, e, d);
	//mpz_mod(tmp, e, tmp);
	//gmp_printf("e*d mod phi(n) = %Zx\n", tmp);

	//////////////////////////////////////Ű�����߰�////////////////////////////////////////////////////
	mpz_sub_ui(tmp, p, 1);
	mpz_mod(dp, d, tmp);

	mpz_sub_ui(tmp, q, 1);
	mpz_mod(dq, d, tmp);

	mpz_invert(qinv, q, p);

	/////////////////////////////////////��ȣȭ////////////////////////////////////////
	mpz_urandomb(m, state, 1024);
	mpz_powm(c, m, e, n);
	gmp_printf("��ȣȭ m = %Zx\n", m);

	/////////////////////////////////////��ȣȭ////////////////////////////////////////
	mpz_powm(m, c, d, n);
	gmp_printf("��ȣȭ m = %Zx\n", m);

	////////////////////////////////////////////////////
	//x, y���ϱ�
	mpz_mod(tmp, c, p);
	mpz_powm(x, tmp, dp, p);
	gmp_printf("x : %Zx\n", x);

	mpz_mod(tmp, c, q);
	mpz_powm(y, tmp, dq, q);
	gmp_printf("y : %Zx\n", y);

	//m ���ϱ�(x-y ���ϱ�)
	mpz_sub(tmp, x, y);
	mpz_mul(tmp, tmp, qinv);
	mpz_mul(tmp, tmp, q);
	mpz_add(tmp, y, tmp);
	mpz_mod(m, tmp, n);
	gmp_printf("�߱��� ������ ���� m : %Zx\n", m);
	
	mpz_clear(p); mpz_clear(q); mpz_clear(n); mpz_clear(e); mpz_clear(d); mpz_clear(tmp); mpz_clear(dp); mpz_clear(dq); mpz_clear(qinv); mpz_clear(m); mpz_clear(c); mpz_clear(x); mpz_clear(y);
	gmp_randclear(state);

	getchar();
	return 0;
}