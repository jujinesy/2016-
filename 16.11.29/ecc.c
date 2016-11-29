#include "ecc.h"




// a mod n = r
int GFP_fast_reduction_p224(mpz_t r,  const mpz_t a, const mpz_t n )
{
	mpz_t prime, s1, s2;
	unsigned long copy_a[14]={0,};
	int i;

	if( (mpz_size(a)<0) || (mpz_size(a)>14) ) 
		return -1;
	
	mpz_init2(prime,224);
	mpz_init2(s1,224);
	mpz_init2(s2,224);

	prime->_mp_d[0] = 1;
	prime->_mp_d[1] = 0;
	prime->_mp_d[2] = 0;
	prime->_mp_d[3] = 0xffffffff;
	prime->_mp_d[4] = 0xffffffff;
	prime->_mp_d[5] = 0xffffffff;
	prime->_mp_d[6] = 0xffffffff;
	prime->_mp_size = 7;

	if( mpz_cmp(prime,n)!=0 ){
		mpz_clear(prime);
		mpz_clear(s1);
		mpz_clear(s2);
	
		return -1;
	}

	for(i=0; i<(int)mpz_size(a); i++) copy_a[i] = a->_mp_d[i];

	//s1
	for(i=0; i<7; i++) s1->_mp_d[i] = copy_a[i];
	s1->_mp_size = 7;

	//s1+s2    , s2 = (c10, c9, c8, c7,0,0,0),
	s2->_mp_d[0] = s2->_mp_d[1] = s2->_mp_d[2] = 0;
	s2->_mp_d[3] = copy_a[7];
	s2->_mp_d[4] = copy_a[8];
	s2->_mp_d[5] = copy_a[9];
	s2->_mp_d[6] = copy_a[10];
	s2->_mp_size = 7;
	mpz_add(s1,s1,s2);

	//s1+s2+s3    s3 = (0, c13, c12, c11,0,0,0)
	s2->_mp_d[0] = s2->_mp_d[1] = s2->_mp_d[2] = s2->_mp_d[6] = 0;
	s2->_mp_d[3] = copy_a[11];
	s2->_mp_d[4] = copy_a[12];
	s2->_mp_d[5] = copy_a[13];
	s2->_mp_size = 7;
	mpz_add(s1,s1,s2);

	//s1+s2+s3-s4   , s4 = (c13, c12, c11, c10, c9, c8, c7),
	s2->_mp_d[0] = copy_a[7];
	s2->_mp_d[1] = copy_a[8];
	s2->_mp_d[2] = copy_a[9];
	s2->_mp_d[3] = copy_a[10];
	s2->_mp_d[4] = copy_a[11];
	s2->_mp_d[5] = copy_a[12];
	s2->_mp_d[6] = copy_a[13];
	s2->_mp_size = 7;
	mpz_sub(s1,s1,s2);

	//s1+s2+s3-s4-s5   , s5 = (0,0,0,0, c13, c12, c11).
	s2->_mp_d[0] = copy_a[11];
	s2->_mp_d[1] = copy_a[12];
	s2->_mp_d[2] = copy_a[13];
	s2->_mp_size = 3;
	mpz_sub(s1,s1,s2);

	if(s1->_mp_size < 0){
		while(s1->_mp_size < 0){
			mpz_add(s1,s1,prime);
		}
	}else{
		while( mpz_cmp(s1,prime)>=0 ){
			mpz_sub(s1,s1,prime);
		}
	}

	mpz_set(r,s1);

	mpz_clear(prime);
	mpz_clear(s1);
	mpz_clear(s2);

	return 0;
}

int GFP_fast_reduction_p256(mpz_t r,  const mpz_t a, const mpz_t n )
{
	return 0;
}




//  r = a+b mod n 
void mpz_add_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n)
{
	mpz_add(r,a,b);
	//if(mpz_cmp(r,n)>=0) mpz_sub(r,r,n);
	mpz_mod(r,r,n);
}
//  r = a-b mod n 
void mpz_sub_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n)
{
	mpz_sub(r,a,b);
	//if(mpz_sgn(r)<0) mpz_add(r,r,n);
	mpz_mod(r,r,n);
}
//  r = a*b mod n 
void mpz_mul_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n)
{
	mpz_mul(r,a,b);
	mpz_mod(r,r,n);
}


void GFP_point_init(GFP_POINT *p)
{
	mpz_init(p->x);
	mpz_init(p->y);
	p->point_at_infinity = 1;
}


void GFP_point_clear(GFP_POINT *p)
{
	mpz_clear(p->x);
	mpz_clear(p->y);
	p->point_at_infinity = 0;
}

void GFP_point_set(GFP_POINT *r,GFP_POINT *p)
{
	if(p->point_at_infinity == 1){
		r->point_at_infinity = 1;
	}else{
		r->point_at_infinity = 0;
		mpz_set(r->x,p->x);
		mpz_set(r->y,p->y);
	}
}


// R=2P
int GFP_affine_doubling(GFP_POINT *r,GFP_POINT *p,const mpz_t coefficient_a,const mpz_t prime)
{
	mpz_t tmp1,tmp2,tmp3,result;

	if(p->point_at_infinity == 1){
		r->point_at_infinity = 1;
		return 0;
	}

	if(p->y->_mp_size == 0){
		return -1;
	}

	mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);mpz_init(result);

	//tmp1 = (3x_1^2+a)*(2y_1)^(-1)
	mpz_mul_mod(tmp1,p->x,p->x,prime);
	mpz_mul_ui(tmp1,tmp1,3);
	mpz_add_mod(tmp1,tmp1,coefficient_a,prime);

	mpz_mul_2exp(tmp2,p->y,1);
	mpz_mod(tmp2,tmp2,prime);
	mpz_invert(tmp2,tmp2,prime);

	mpz_mul_mod(tmp1,tmp1,tmp2,prime);

	// x3
	mpz_mul(tmp2,tmp1,tmp1);
	mpz_mod(tmp2,tmp2,prime);
	mpz_mul_2exp(tmp3,p->x,1);
	mpz_mod(tmp3,tmp3,prime);

	mpz_sub_mod(result,tmp2,tmp3,prime);


	// y3
	mpz_sub(tmp2,p->x,result);
	mpz_mul_mod(tmp3,tmp1,tmp2,prime);
	mpz_sub_mod(r->y,tmp3,p->y,prime);
	mpz_set(r->x,result);


	mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);mpz_clear(result);


	return 0;
}
// R=P+Q
int GFP_affine_addition(GFP_POINT *r,GFP_POINT *p,GFP_POINT *q,const mpz_t coefficient_a,const mpz_t prime)
{
	mpz_t tmp1,tmp2,tmp3,result;

	if(p->point_at_infinity == 1){

		GFP_point_set(r,q);
		return 0;
	}
	
	if(q->point_at_infinity == 1){

		GFP_point_set(r,p);
		return 0;
	}

	if( mpz_cmp(p->x,q->x)==0 ){

		if( mpz_cmp(p->y,q->y)==0 ){
			return GFP_affine_doubling(r,p,coefficient_a,prime);			
		}

		r->point_at_infinity = 1;

		return 0;

	}else{
		mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);mpz_init(result);

		//±â¿ï±â
		mpz_sub_mod(tmp1,q->x,p->x,prime);
		mpz_invert(tmp1,tmp1,prime);
		mpz_sub(tmp2,q->y,p->y);
		mpz_mul_mod(tmp1,tmp1,tmp2,prime);

		//x3
		mpz_mul_mod(tmp2,tmp1,tmp1,prime);
		mpz_sub(tmp2,tmp2,p->x);
		mpz_sub_mod(result,tmp2,q->x,prime); 

		//y3
		mpz_sub(tmp2,p->x,result);
		mpz_mul_mod(tmp2,tmp1,tmp2,prime);
		mpz_sub_mod(r->y,tmp2,p->y,prime);
		mpz_set(r->x,result);

		mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);mpz_clear(result);
	}

	return 0;
}


// R=kP LtoR Binary
int GFP_LtoR_binary(GFP_POINT *r,mpz_t k,GFP_POINT *p,const mpz_t coefficient_a,const mpz_t prime)
{
	int i,j;
	GFP_POINT result;

	if( k->_mp_size<0 ) 
		return -1;

	GFP_point_init(&result);
	
	for(i=mpz_size(k)-1 ; i>=0 ; i--)
	{
		for(j=31 ; j>=0 ; j--)
		{
			GFP_affine_doubling(&result,&result,coefficient_a,prime);	

			if( k->_mp_d[i] & (1<<j) )
			{
				GFP_affine_addition(&result,&result,p,coefficient_a,prime);
			}
		}
	}

	GFP_point_set(r,&result);

	GFP_point_clear(&result);

	return 0;
}


// R=kP LtoR NAF
int GFP_naf_recording(NAF_RECORDING *rk, mpz_t k)
{
	int i=0;
	unsigned long tmp;

	if( k->_mp_size<0 ) 
		return -1;

	while( mpz_cmp_ui(k,0)!=0 )
	{
		tmp = k->_mp_d[0];
		if( tmp&1 ){
			rk->naf_scalar[i] = 2 - (tmp & 0x3);

			if(rk->naf_scalar[i]>0){
				mpz_sub_ui(k,k,rk->naf_scalar[i]);
			}else{
				mpz_add_ui(k,k,(-1)*rk->naf_scalar[i]);
			}
		}else{
			rk->naf_scalar[i] = 0;
		}
		mpz_tdiv_q_2exp(k,k,1);
		i++;
	}

	rk->naf_len = i;

	return 0;
}

int GFP_LtoR_NAF(GFP_POINT *r,NAF_RECORDING *k,GFP_POINT *p,const mpz_t coefficient_a,const mpz_t prime)
{
	int i;
	GFP_POINT result,neg_p;

	if( k->naf_len>1024 ) 
		return -1;

	GFP_point_init(&result);
	GFP_point_init(&neg_p);
	GFP_point_set(&neg_p,p);
	mpz_sub(neg_p.y,prime,neg_p.y);
	//neg_p.y->_mp_size = neg_p.y->_mp_size*(-1);
	
	for(i=k->naf_len-1 ; i>=0 ; i--)
	{
		GFP_affine_doubling(&result,&result,coefficient_a,prime);	

		if( k->naf_scalar[i] != 0 )
		{
			if( k->naf_scalar[i] > 0 )
				GFP_affine_addition(&result,&result,p,coefficient_a,prime);
			else
				GFP_affine_addition(&result,&result,&neg_p,coefficient_a,prime);  
		}		
	}

	GFP_point_set(r,&result);

	GFP_point_clear(&result);

	return 0;
}