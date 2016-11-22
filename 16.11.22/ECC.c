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






// r = a + b mod n
void mpz_add_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n) {
	mpz_add(r, a, b);
	if (mpz_cmp(r, n) >= 0)  mpz_sub(r, r, n);
}

// r = a - b mod n
void mpz_sub_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n) {
	mpz_sub(r, a, b);
	if (mpz_sgn(r) < 0)  mpz_add(r, r, n);
}

// r = a * b mod n
void mpz_mul_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n) {
	mpz_mul(r, a, b);
	mpz_mod(r, r, n);
}



////cr  ap
//int GFP_affine_doubling(GFP_POINT *r, GFP_POINT *p, const mpz_t coefficient_a, const mpz_t prime)
//{
//	mpz_t t1, t2;
//
//	mpz_init(t1);
//	mpz_init(t2);
//
//	if ((mpz_cmp_ui(p->x, 0) == 0) && (mpz_cmp(p->y, 0) == 0)) {
//		mpz_set_ui(r->x, 0);
//		mpz_set_ui(r->y, 0);
//		r->point_at_infinity = 1;
//		return;
//	}
//
//	if (mpz_cmp_ui(p->y, 0) == 0) {
//		//  error
//		return;
//	}
//
//	mpz_mul_mod(t1, p->x, p->x, prime);
//	mpz_set_ui(r->x, 3);
//	mpz_mul_mod(t1, t1, r->x, prime);   //  3x1^2
//	mpz_add_mod(t1, t1, coefficient_a, prime);
//	mpz_set_ui(r->x, 2);
//	mpz_mul_mod(r->x, r->x, p->y, prime);
//	mpz_invert(r->x, r->x, prime);
//	mpz_mul_mod(t1, t1, r->x, prime);
//
//	mpz_mul_mod(t2, t1, t1, prime);
//
//	mpz_set_ui(r->x, 2);
//	mpz_mul_mod(r->x, r->x, p->x, prime);
//	mpz_sub_mod(r->x, t2, r->x, prime);
//
//	mpz_sub_mod(r->y, p->x, r->x, prime);
//	mpz_mul_mod(r->y, r->y, t1, prime);
//	mpz_sub_mod(r->y, r->y, p->y, prime);
//
//	mpz_clear(t1);
//	mpz_clear(t2);
//}
//R = 2P
int GFP_affine_doubling(GFP_POINT *r, GFP_POINT *p, const mpz_t coeffcient_a, const mpz_t prime) {
	mpz_t tmp, tmp2, t1, t2, tmpx, tmpy;
	mpz_init(tmp);

	if (p->point_at_infinity == 1) {
		r->point_at_infinity = 1;
		mpz_clear(tmp);
		return 0;
	}
	else if (mpz_cmp(p->y, tmp) == 0) {
		mpz_clear(tmp);
		return -1;
	}
	mpz_init(tmp2);
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(tmpx);
	mpz_init(tmpy);

	mpz_set(tmpx, p->x);
	mpz_set(tmpy, p->y);
	//t1 = x mod prime
	mpz_mod(t1, p->x, prime);

	//t1 = 3x^2 + a
	mpz_mul(t1, t1, t1);
	mpz_mul_ui(t1, t1, 3);
	mpz_add(t1, t1, coeffcient_a);

	//t2 = 2y^-1
	mpz_mul_ui(t2, p->y, 2);
	mpz_invert(t2, t2, prime);

	//tmp = (3x^2+a)/2y
	//tmp2 = {(3x^2+a)/2y}^2
	//t1 = 2x
	mpz_mul(tmp, t1, t2);
	mpz_mul(tmp2, tmp, tmp);
	mpz_mul_ui(t1, p->x, 2);
	mpz_sub(tmp2, tmp2, t1);
	mpz_mod(r->x, tmp2, prime);

	////////////////////////////
	//t1 = x1 - x3
	mpz_sub(t1, tmpx, r->x);
	//
	mpz_mul(tmp, tmp, t1);
	mpz_sub(tmp2, tmp, p->y);
	mpz_mod(r->y, tmp2, prime);

	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(tmpx);
	mpz_clear(tmpy);
	return 0;
}

////cr  ap bq
//int GFP_affine_addition(GFP_POINT *r, GFP_POINT *p, GFP_POINT *q, const mpz_t coefficient_a, const mpz_t prime)
//{
//	mpz_t t;
//	mpz_init(t);
//
//	if ((mpz_cmp_ui(p->x, 0) == 0) && (mpz_cmp(p->y, 0) == 0)) {
//		mpz_set(r->x, q->x);
//		mpz_set(r->y, q->x);
//		r->point_at_infinity = q->point_at_infinity;
//		return;
//	}
//
//	if ((mpz_cmp_ui(q->x, 0) == 0) && (mpz_cmp(q->y, 0) == 0)) {
//		mpz_set(r->x, p->x);
//		mpz_set(r->y, p->x);
//		r->point_at_infinity = p->point_at_infinity;
//		return;
//	}
//
//	if (mpz_cmp(p->x, q->x) == 0) {
//		if (mpz_cmp(p->y, q->y) == 0) {
//			GFP_affine_doubling(r, p, coefficient_a, prime);
//			return;
//		}
//		else {
//			mpz_set(r->x, q->x);
//			mpz_set(r->y, q->x);
//			r->point_at_infinity = q->point_at_infinity;
//			return;
//		}
//	}
//
//	mpz_sub_mod(r->x, q->y, p->y, prime);
//	mpz_sub_mod(t, q->x, p->x, prime);
//	mpz_invert(t, t, prime);
//
//	mpz_mul_mod(t, r->x, t, prime);
//	mpz_mul_mod(r->x, t, t, prime);
//
//	mpz_sub_mod(r->x, r->x, p->x, prime);
//	mpz_sub_mod(r->x, r->x, q->x, prime);
//
//	//
//	mpz_sub_mod(r->y, p->x, r->x, prime);
//	mpz_mul_mod(r->y, t, r->y, prime);
//	mpz_sub_mod(r->y, r->y, p->y, prime);
//
//	mpz_clear(t);
//}
//R = P + Q
int GFP_affine_addition(GFP_POINT *r, GFP_POINT *p, GFP_POINT *q, const mpz_t coeffcient_a, const mpz_t prime) {

	mpz_t tmp, tmp2, t1, t2, tmpx, tmpy;
	mpz_init(tmp);

	if (p->point_at_infinity == 1) {
		mpz_set(r->x, q->x);
		mpz_set(r->y, q->y);
		r->point_at_infinity = q->point_at_infinity;
		mpz_clear(tmp);
		return 0;
	}
	else if (q->point_at_infinity == 1) {
		mpz_set(r->x, p->x);
		mpz_set(r->y, p->y);
		r->point_at_infinity = p->point_at_infinity;
		mpz_clear(tmp);
		return 0;
	}

	if (mpz_cmp(p->x, q->x) == 0) {
		if (mpz_cmp(p->y, q->y) == 0) {
			r->point_at_infinity = 1;
			mpz_clear(tmp);
			return 0;
		}
		else {
			GFP_affine_doubling(r, p, coeffcient_a, prime);
			mpz_clear(tmp);
			return 0;
		}
	}

	mpz_init(tmp2);
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(tmpx);
	mpz_init(tmpy);

	mpz_set(tmpx, p->x);
	mpz_set(tmpy, p->y);

	// t1 = y2 - y1
	mpz_sub(t1, q->y, p->y);
	// t2 = (x2 - x1)^-1
	mpz_sub(t2, q->x, p->x);

	mpz_mod(t2, t2, prime);
	mpz_invert(t2, t2, prime);
	// tmp = (y2 - y1) * (x2 - x1)^-1
	mpz_mul(tmp, t1, t2);
	mpz_mod(tmp, tmp, prime);

	// tmp2 = {(y2 - y1) / (x2 - x1)}^2
	mpz_mul(tmp2, tmp, tmp);
	mpz_mod(tmp2, tmp2, prime);

	// tmp2 = {(y2 - y1) / (x2 - x1)}^2 -x1
	mpz_sub(tmp2, tmp2, p->x);
	// tmp2 = {(y2 - y1) / (x2 - x1)}^2 -x1 -x2
	mpz_sub(tmp2, tmp2, q->x);
	mpz_mod(r->x, tmp2, prime);

	//////////////////////////////
	//t1 = x1 - x3
	mpz_sub(t1, tmpx, r->x);
	//
	mpz_mul(tmp2, tmp, t1);
	mpz_sub(tmp2, tmp2, p->y);
	mpz_mod(r->y, tmp2, prime);

	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(tmpx);
	mpz_clear(tmpy);
	return 0;
}

int GFP_point_init(GFP_POINT *p)
{
	mpz_init(p->x);
	mpz_init(p->y);
	p->point_at_infinity = 1;
}


int GFP_point_clear(GFP_POINT *p)
{
	mpz_clear(p->x);
	mpz_clear(p->y);
	p->point_at_infinity = 0;
}


void GFP_mod(mpz_t c, mpz_t a, mpz_t n)
//  c= a mod n 
{
#ifndef _224_FAST_REDUCTION_
	mpz_mod(c, a, n);
#else
	GFP_P224_fast_reduction(c, a, n);
#endif

}

void GFP_add_mod(mpz_t c, mpz_t a, mpz_t b, mpz_t n)
{
	mpz_add(c, a, b);
	GFP_mod(c, c, n);

}

void GFP_sub_mod(mpz_t c, mpz_t a, mpz_t b, mpz_t n)
// c=a-b mod n
{
	mpz_sub(c, a, b);
	GFP_mod(c, c, n);
}

void GFP_mul_mod(mpz_t c, mpz_t a, mpz_t b, mpz_t n)
// c=a*b mod n
{
	mpz_mul(c, a, b);
	GFP_mod(c, c, n);
}

void GFP_invert(mpz_t c, mpz_t a, mpz_t n)
// c=a^(-1) mod n
{
	mpz_invert(c, a, n);

}



//int GFP_LtoR_binary(GFP_POINT *r, mpz_t k, GFP_POINT *p, const mpz_t coefficient_a, const mpz_t prime)
//{
//	int i, j;
//	GFP_POINT *result, *tmp;
//	GFP_point_init(&result);
//	GFP_point_init(&tmp);
//
//	if (k->_mp_size < 0)
//		return -1;
//
//	GFP_point_init(&result); GFP_point_init(&tmp);
//
//	for (i = mpz_size(k) - 1; i <= 0; i--)
//	{
//		for (j = 31; j >= 0; j--)
//		{
//			GFP_affine_doubling(&tmp, &result, coefficient_a, prime);
//
//
//			mpz_set(result->x, tmp->x);
//			mpz_set(result->y, tmp->y);
//			tmp->point_at_infinity = result->point_at_infinity;
//
//
//			if (k->_mp_d[i] & (1 << j))
//			{
//				GFP_affine_addition(&tmp, &result, p, coefficient_a, prime);
//
//				mpz_set(result->x, tmp->x);
//				mpz_set(result->y, tmp->y);
//				tmp->point_at_infinity = result->point_at_infinity;
//
//			}
//		}
//	}
//
//	//GFP_point_set(r, &result);
//	mpz_set(r->x, result->x);
//	mpz_set(r->y, result->y);
//	result->point_at_infinity = r->point_at_infinity;
//
//
//	mpz_clear(result); mpz_clear(tmp);
//
//	return 0;
//
//}
//R = kP LtoR binary
int GFP_LtoR_binary(GFP_POINT *r, mpz_t k, GFP_POINT *p, const mpz_t coeffcient_a, const mpz_t prime) {

	int i, j;
	GFP_POINT tmp;

	if (k->_mp_size < 0) return -1;

	GFP_point_init(&tmp);

	for (i = mpz_size(k) - 1; i >= 0; i--) {
		for (j = 31; j >= 0; j--) {
			GFP_affine_doubling(&tmp, &tmp, coeffcient_a, prime);
			if (k->_mp_d[i] & (1 << j)) {
				GFP_affine_addition(&tmp, &tmp, p, coeffcient_a, prime);
			}
		}
	}
	mpz_set(r->x, tmp.x);
	mpz_set(r->y, tmp.y);
	r->point_at_infinity = tmp.point_at_infinity;

	GFP_point_clear(&tmp);
	return 0;
}

int GFP_naf_recording(NAF_RECORDING *rk, mpz_t k) //3.30
{
	int i = 0;
	unsigned long tmp;

	while (mpz_cmp_ui(k, 0) != 0)
	{
		tmp = k->_mp_d[0];
		if (tmp & 1) {
			rk->naf_scalar[i] = 2 - (tmp & 0x3);
			if (rk->naf_scalar[i]>0) {
				mpz_sub_ui(k, k, rk->naf_scalar[i]);
			}
			else {
				mpz_add_ui(k, k, (-1)*rk->naf_scalar[i]);
			}
		}
		else {
			rk->naf_scalar[i] = 0;
		}
		mpz_tdiv_q_2exp(k, k, 1);
		i = i + 1;
	}
	rk->naf_len = i;
}


int GFP_LtoR_NAR(GFP_POINT *r, NAF_RECORDING *k, GFP_POINT *p, const mpz_t coefficient_a, const mpz_t prime) {
	int i, j;
	GFP_POINT result, tmp, neg_p;

	GFP_point_init(&result);
	GFP_point_init(&tmp);
	GFP_point_init(&neg_p);

	mpz_set(neg_p.x, p->x);
	mpz_set(neg_p.y, p->y);
	neg_p.point_at_infinity = p->point_at_infinity;

	mpz_sub(neg_p.y, prime, neg_p.y);

	for (i = k->naf_len - 1; i >= 0; i--) {
		GFP_affine_doubling(&tmp, &result, coefficient_a, prime);
		mpz_set(result.x, tmp.x);
		mpz_set(result.y, tmp.y);
		result.point_at_infinity = tmp.point_at_infinity;

		if (k->naf_scalar[i] != 0) {
			if (k->naf_scalar[i] > 0) {
				GFP_affine_addition(&tmp, &result, p, coefficient_a, prime);
				mpz_set(result.x, tmp.x);
				mpz_set(result.y, tmp.y);
				result.point_at_infinity = tmp.point_at_infinity;
			}
			else {
				GFP_affine_addition(&tmp, &result, &neg_p, coefficient_a, prime);
				mpz_set(result.x, tmp.x);
				mpz_set(result.y, tmp.y);
				result.point_at_infinity = tmp.point_at_infinity;
			}
		}
	}

	GFP_point_clear(&result); GFP_point_clear(&tmp); GFP_point_clear(&neg_p);
	return 0;
}