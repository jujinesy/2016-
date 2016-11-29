#pragma comment(lib, "../include/gmp.lib") 
#include "../include/gmp.h"
#include "sha2.h"


//#define HASH_ALGO_SHA1     0
#define HASH_ALGO_SHA224   1
#define HASH_ALGO_SHA256   2
#define HASH_ALGO_SHA384   3
#define HASH_ALGO_SHA512   4


typedef struct _ECC_HASH_ALGO{
	void	(*hash)(const unsigned char *M, unsigned int len, unsigned char *Hash);
	int  	digest_byte_len;
	int 	block_byte_len;
}ECC_HASH_ALGO;

static const ECC_HASH_ALGO ecc_hash_algo[10]={
	{0,0,0},
	{sha224,SHA224_DIGEST_SIZE,SHA224_BLOCK_SIZE},
	{sha256,SHA256_DIGEST_SIZE,SHA256_BLOCK_SIZE},
	{sha384,SHA384_DIGEST_SIZE,SHA384_BLOCK_SIZE},
	{sha512,SHA512_DIGEST_SIZE,SHA512_BLOCK_SIZE},
};

typedef struct _GFP_POINT_
{
	mpz_t x;
	mpz_t y;
	int   point_at_infinity;
}GFP_POINT; 



typedef struct _NAF_RECORDING_
{
	char naf_scalar[1024];
	int  naf_len;
}NAF_RECORDING; 




typedef struct _ECC_DOMAIN_PARAM_{
	int       curve_num;
	mpz_t     prime;
	mpz_t     coeff_a;
	mpz_t     coeff_b;
	GFP_POINT base;
	mpz_t     order;
	mpz_t     cofactor;
}ECC_DOMAIN_PARAM;




// a mod n = r
int GFP_fast_reduction_p224(mpz_t r,  const mpz_t a, const mpz_t n );
int GFP_fast_reduction_p256(mpz_t r,  const mpz_t a, const mpz_t n );

//  r = a+b mod n 
void mpz_add_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n);
//  r = a-b mod n 
void mpz_sub_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n);
//  r = a*b mod n 
void mpz_mul_mod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n);



// ECC computation
void GFP_point_init(GFP_POINT *p);
void GFP_point_clear(GFP_POINT *p);
void GFP_point_set(GFP_POINT *r,GFP_POINT *p);

// R=2P
int GFP_affine_doubling(GFP_POINT *r,GFP_POINT *p,const mpz_t coefficient_a,const mpz_t prime); 
// R=P+Q
int GFP_affine_addition(GFP_POINT *r,GFP_POINT *p,GFP_POINT *q,const mpz_t coefficient_a,const mpz_t prime); 



// R=kP LtoR binary
int GFP_LtoR_binary(GFP_POINT *r,mpz_t k,GFP_POINT *p,const mpz_t coefficient_a,const mpz_t prime); 

// R=kP LtoR NAF
int GFP_naf_recording(NAF_RECORDING *rk, mpz_t k);
int GFP_LtoR_NAF(GFP_POINT *r,NAF_RECORDING *k,GFP_POINT *p,const mpz_t coefficient_a,const mpz_t prime); 



//ECDSA keypair gen
int ECDSA_key_pair_gen(mpz_t prikey,//비밀키
	GFP_POINT *pubkey, //공개키
	GFP_POINT *p, //생성원
	mpz_t order, //타원곡선 위수
	const mpz_t coefficient_a, //타원곡선계수
	const mpz_t prime);//유한체 소수

int ECDSA_sign_gen(mpz_t r, //서명값 r
	mpz_t s, //서명값 s
	unsigned char *message, //메시지
	int message_len, //메시지 길이
	mpz_t prikey, //비밀키
	mpz_t order, //타원곡선 위수
	int hash_algo_info, //해시함수 정보
	GFP_POINT *p, //생성원
	const mpz_t coefficient_a, //타원곡선계수
	const mpz_t prime);//유한체 소수

int ECDSA_sign_ver(mpz_t r, //서명값 r
	mpz_t s, //서명값 s
	unsigned char *message, //메시지
	int message_len, //메시지 길이
	GFP_POINT *pubkey, //공개키
	mpz_t order, //타원곡선 위수
	int hash_algo_info, //해시함수 정보
	GFP_POINT *p, //생성원
	const mpz_t coefficient_a, //타원곡선계수
	const mpz_t prime);//유한체 소수

int ostr2mpz(mpz_t a, const unsigned char *ostr, const int ostrlen);