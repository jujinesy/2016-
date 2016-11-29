#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ecc.h"



clock_t elapsed;
float sec;

#define START_WATCH \
{\
 elapsed = -clock(); \
}\
 
#define STOP_WATCH \
{\
 elapsed += clock();\
 sec = (float)elapsed/CLOCKS_PER_SEC;\
}\
 
#define PRINT_TIME(qstr) \
{\
 printf("\n[%s: %.5f s]\n",qstr,sec);\
}\



void test_div_mul()
{
	mpz_t a,b,c,d,q,r;
	gmp_randstate_t state;
	int i;

	mpz_init(a);mpz_init(b);mpz_init(c);mpz_init(d);
	mpz_init(q);mpz_init(r);
	gmp_randinit_default(state);

	i = mpz_size(a); 


	//a »ý¼º
	mpz_urandomb(a,state,1024);
	mpz_urandomb(b,state,1024);
	mpz_urandomb(d,state,1024);


	START_WATCH;
	for(i=0;i<100000;i++){
		mpz_urandomb(a,state,1024);
		mpz_urandomb(b,state,1024);
		mpz_urandomb(d,state,1024);

	}
	STOP_WATCH;
	PRINT_TIME("base time");


	START_WATCH;
	for(i=0;i<100000;i++){
		mpz_add(c,a,b);

	}
	STOP_WATCH;
	PRINT_TIME("1024bit add");

	START_WATCH;
	for(i=0;i<100000;i++){
		mpz_mul(c,a,b);

	}
	STOP_WATCH;
	PRINT_TIME("1024bit mul");


	START_WATCH;
	for(i=0;i<100000;i++){
		mpz_mod(r,c,d);
	}
	STOP_WATCH;
	PRINT_TIME("1024bit div");
	


	mpz_clear(a);mpz_clear(b);mpz_clear(c);mpz_clear(d);
	mpz_clear(q);mpz_clear(r);
	gmp_randclear(state);

}

void test_p224_mod()
{
	mpz_t a,b,c,d,e,prime;
	gmp_randstate_t state;
	int i;

	mpz_init(a);mpz_init(b);mpz_init(c);mpz_init(d);mpz_init(e);
	mpz_init2(prime,224);
	gmp_randinit_default(state);

	prime->_mp_d[0] = 1;
	prime->_mp_d[1] = 0;
	prime->_mp_d[2] = 0;
	prime->_mp_d[3] = 0xffffffff;
	prime->_mp_d[4] = 0xffffffff;
	prime->_mp_d[5] = 0xffffffff;
	prime->_mp_d[6] = 0xffffffff;
	prime->_mp_size = 7;

	mpz_urandomm(a,state,prime);
	mpz_urandomm(b,state,prime);

	mpz_mul(c,a,b);

	START_WATCH;
	for(i=0;i<100000;i++){
		mpz_mod(d,c,prime);

	}
	STOP_WATCH;
	PRINT_TIME("gen mod");

	START_WATCH;
	for(i=0;i<100000;i++){
		GFP_fast_reduction_p224(e,c,prime);

	}
	STOP_WATCH;
	PRINT_TIME("fast mod");	

	

	if(mpz_cmp(d,e)==0){
		printf("true");
	}
	else{
		printf("false");
	}

	
}


void ecc_affine_add_doubling_test()
{

	GFP_POINT p,q,r,pq,p2;
	mpz_t prime,co_a;

	GFP_point_init(&p);
	GFP_point_init(&q);
	GFP_point_init(&r);
	GFP_point_init(&pq);
	GFP_point_init(&p2);
	mpz_init(prime);
	mpz_init(co_a);



	//coefficeit_a
	mpz_set_str(co_a, "fffffffffffffffffffffffffffffffefffffffffffffffffffffffe",16);

	//prime
	mpz_set_str(prime,"ffffffffffffffffffffffffffffffff000000000000000000000001",16);


	//p 
	mpz_set_str(p.x,"b70e0cbd6bb4bf7f321390b94a03c1d356c21122343280d6115c1d21",16);
	mpz_set_str(p.y,"bd376388b5f723fb4c22dfe6cd4375a05a07476444d5819985007e34",16);
	p.point_at_infinity=0;

	//q
	mpz_set_str(q.x,"3e54b1067d79dab2bd159a61f03e7bcb62851c5d3418a8bc5cad16aa",16);
	mpz_set_str(q.y,"5986453b559d4ce6b24a4f10a87efe19eaf49a6db877bda6b791adb3",16);
	q.point_at_infinity=0;


	//addition
	GFP_affine_addition(&r,&p,&q,co_a,prime); 


	//doubling

	GFP_affine_doubling(&p2,&p,co_a,prime);

	gmp_printf("add rx = %Zx\n",r.x);
	gmp_printf("add ry = %Zx\n",r.y);

	gmp_printf("dbl rx = %Zx\n",p2.x);
	gmp_printf("dbl ry = %Zx\n",p2.y);


	GFP_point_clear(&p);
	GFP_point_clear(&q);
	GFP_point_clear(&r);
	GFP_point_clear(&pq);
	GFP_point_clear(&p2);
	mpz_clear(prime);
	mpz_clear(co_a);
}

void test_ecc_scalar_mul_test()
{
	
	GFP_POINT p,q,r;
	mpz_t prime,co_a,k;
	NAF_RECORDING rk;

	GFP_point_init(&p);
	GFP_point_init(&q);
	GFP_point_init(&r);
	mpz_init(prime);
	mpz_init(k);
	mpz_init(co_a);


	//coefficeit_a
	mpz_set_str(co_a, "fffffffffffffffffffffffffffffffefffffffffffffffffffffffe",16);

	//prime
	mpz_set_str(prime,"ffffffffffffffffffffffffffffffff000000000000000000000001",16);


	//p 
	mpz_set_str(p.x,"b70e0cbd6bb4bf7f321390b94a03c1d356c21122343280d6115c1d21",16);
	mpz_set_str(p.y,"bd376388b5f723fb4c22dfe6cd4375a05a07476444d5819985007e34",16);
	p.point_at_infinity=0;

	//k
	mpz_set_str(k,"3354543463756857898797965432",16);



	GFP_LtoR_binary(&q,k,&p,co_a,prime); 

    GFP_naf_recording(&rk,k);
	GFP_LtoR_NAF(&r,&rk,&p,co_a,prime); 


	gmp_printf("add rx = %Zx\n",q.x);
	gmp_printf("add ry = %Zx\n",q.y);

	gmp_printf("dbl rx = %Zx\n",r.x);
	gmp_printf("dbl ry = %Zx\n",r.y);


	GFP_point_clear(&p);
	GFP_point_clear(&q);
	GFP_point_clear(&r);
	mpz_clear(prime);
	mpz_clear(k);
	mpz_clear(co_a);

}


void main()
{
	//test_div_mul();
	//test_p224_mod();
	//ecc_affine_add_doubling_test();
	test_ecc_scalar_mul_test();
}

