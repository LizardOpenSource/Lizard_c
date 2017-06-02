#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <immintrin.h>
#include <stdint.h>

// compile and run example:
// $ gcc Lizard.c -o Lizard -O2
// $ ./Lizard

//////////////////////////////////  PARAMETER CHOICES  /////////////////////////////////////
//#define PARAMS_Classical								// classical parameter
//#define PARAMS_Recommended							// recommended parameter
//#define PARAMS_Classical_Plaintext_32bit				// classical parameter with 32-bit plaintext
#define PARAMS_CCA									// CCA parameter
//#define PARAMS_RING									// Ring version parameter

//////////////////////////////////  PARAMETER SETS  ////////////////////////////////////////

#define iter 100 		// iteration number for keygen & EncDec test
#define testnum 1000	// repeatetion number of Enc Dec procedure in a single iteration

#ifdef PARAMS_Classical
#define PARAMNAME "Classical"
#define LWE_N 544		// LWE dim
#define LWE_M 940		// LWR dim (Number of LWE samples in pk)
#define LWE_L 256
#define LOG_Q 10
#define _16_LOG_Q 6
#define LOG_P 8
#define RD_ADD 0x80 	// 2^(15 - LOG_P)
#define RD_AND 0xff00
#define LOG_T 1
#define _16_LOG_T 15
#define T 2
#define DEC_ADD 0x4000	// 2^(15 - LOG_T)
#define NOISE_D1		// standard deviation for discrete gaussian distribution
#define HR 128			// Hamming weight of coefficient vector r
#endif

#ifdef PARAMS_Recommended
#define PARAMNAME "Recommended"
#define LWE_N 608		// LWE dim
#define LWE_M 960		// LWR dim (Number of LWE samples in pk)
#define LWE_L 256
#define LOG_Q 10
#define _16_LOG_Q 6
#define LOG_P 8
#define RD_ADD 0x80 	// 2^(15 - LOG_P)
#define RD_AND 0xff00
#define LOG_T 1
#define _16_LOG_T 15
#define T 2
#define DEC_ADD 0x4000	// 2^(15 - LOG_T)
#define NOISE_D2		// standard deviation for discrete gaussian distribution
#define HR 128			// Hamming weight of coefficient vector r
#endif

#ifdef PARAMS_Classical_Plaintext_32bit
#define PARAMNAME "Classical_32bit"
#define LWE_N 544		// LWE dim
#define LWE_M 840		// LWR dim (Number of LWE samples in pk)
#define LWE_L 32
#define LOG_Q 10
#define _16_LOG_Q 6
#define LOG_P 8
#define RD_ADD 0x80 	// 2^(15 - LOG_P)
#define RD_AND 0xff00
#define LOG_T 1
#define _16_LOG_T 15
#define T 2
#define DEC_ADD 0x4000
#define NOISE_D1		// standard deviation for discrete gaussian distribution
#define HR 128			// Hamming weight of coefficient vector r
#endif

#ifdef PARAMS_CCA
#define PARAMNAME "CCA"
#define LWE_N 608		// LWE dim
#define LWE_M 1024		// LWR dim (Number of LWE samples in pk)
#define LWE_L 384
#define LOG_Q 10
#define _16_LOG_Q 6
#define LOG_P 8
#define RD_ADD 0x80 	// 2^(15 - LOG_P)
#define RD_AND 0xff00
#define HR 128
#define LOG_T 1
#define _16_LOG_T 15
#define T 2
#define DEC_ADD 0x4000
#define NOISE_D2		// standard deviation for discrete gaussian distribution
#endif

#ifdef PARAMS_RING
#define PARAMNAME "RingLizard"
#define LWE_N 1024		// LWE dim
#define LWE_M 1024		// LWR dim (Number of LWE samples in pk)
#define LWE_L 1024
#define LOG_Q 10
#define _16_LOG_Q 6
#define LOG_P 8
#define RD_ADD 0x80 	// 2^(15 - LOG_P)
#define RD_AND 0xff00
#define HR 128
#define LOG_T 1
#define _16_LOG_T 15
#define T 2
#define DEC_ADD 0x4000
#define NOISE_D2		// standard deviation for discrete gaussian distribution
#endif

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// NOISE DISTRIBUTION /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

// error distribution for LWE instances in public key generation:
// store probability table (CDF_TABLE)
// sample a single error from (RANDBITS)-length of random binary string

//A discrete error distribution close to the discrete Gaussian distribution with sigma = 5.9883 w.r.t. Renyi div = 1.00694
#ifdef NOISE_D1
#define SAMPLE_DG Sample_D1
const uint16_t CDF_TABLE[8] = {42, 120, 180, 219, 240, 250, 254, 255}; // out of [0, 255]
const size_t RANDBITS = 9;
const size_t TABLE_LENGTH = 8;

uint16_t Sample_D1(){
	uint16_t rnd = rand() & 0x00ff;
	uint16_t sign = rand() & 0x01;
	uint16_t sample = 0;
	for(size_t i = 0; i < TABLE_LENGTH - 1; ++i){
		sample += (CDF_TABLE[i] - rnd) >> 15;
	}
	sample = ((-sign) ^ sample) + sign;
	return sample;
}
#endif

//A discrete error distribution close to the discrete Gaussian distribution with sigma = 5.6264 w.r.t. Renyi div = 1.00302
#ifdef NOISE_D2
#define SAMPLE_DG Sample_D2
const uint16_t CDF_TABLE[8] = {90, 255, 377, 452, 489, 505, 510, 511}; // out of [0, 511]
const size_t RANDBITS = 10;
const size_t TABLE_LENGTH = 8;

uint16_t Sample_D2(){
	uint16_t rnd = rand() & 0x01ff;
	uint16_t sign = rand() & 0x01;
	uint16_t sample = 0;
	for(size_t i = 0; i < TABLE_LENGTH - 1; ++i){
		sample += (CDF_TABLE[i] - rnd) >> 15;
	}
	sample = ((-sign) ^ sample) + sign;
	return sample;
}
#endif

clock_t start, finish, elapsed1, elapsed2;

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// CPA SCHEME //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

// Every variables are defined over uint16_t
// For smaller modulus operations, fill the zeros in LSBs of uint16_t
// e.g. when q=2^10, pad 6 bits in LSBs during pk generation

uint16_t sk_CPA[LWE_N * LWE_L];

struct {
	uint16_t A[LWE_M * LWE_N];
	uint16_t B[LWE_M * LWE_L];
} pk_CPA;

uint16_t ptx_CPA[LWE_L];
uint16_t backup_CPA[LWE_L];

struct { // n-dimensional LWE ctxtext modulo q
	uint16_t a[LWE_N];		// array of n integers modulo q
	uint16_t b[LWE_L];		// integer modulo q
} ctx_CPA;

// generate a random (LWE_M * LWE_N) matrix A of public key
void gen_A_CPA(){
	for(int i = 0; i < LWE_M; ++i){
		uint16_t* pk_Ai = pk_CPA.A + LWE_N * i;
		for(int j = 0; j < LWE_N; ++j){
			pk_Ai[j] = rand() << _16_LOG_Q;
		}
	}
}

// generate a (LWE_M * LWE_L) LWE error matrix E in generation of public key
void gen_E_CPA(){
	for(int i = 0; i < LWE_M; ++i){
		uint16_t* pk_Bi = pk_CPA.B + LWE_L * i;
		for(int j = 0; j < LWE_L; ++j){
			pk_Bi[j] = SAMPLE_DG() << _16_LOG_Q;
		}
	}
}

// generate secret key of CPA & CCA version Lizard
void gen_sk_CPA(){
	for(int i = 0; i < LWE_L; ++i){
		uint16_t* sk_i = sk_CPA + LWE_N * i;
		for(int j = 0; j < LWE_N; ++j){
			sk_i[j] = (rand() & 0x01) + (rand() & 0x01) - 1;
		}
	}
}

// generate secret key and public key of CPA & CCA version Lizard
// Both CPA and CCA Lizard schemes use CPA key generation: output sk_CPA and pk_CPA
void Keygen_CPA(){
    elapsed1 = clock();
    for(int l = 0; l < iter; ++l){
	gen_sk_CPA();
	gen_A_CPA();
	gen_E_CPA();

	for(int i = 0; i < LWE_M; ++i){
		uint16_t* A_i = pk_CPA.A + LWE_N * i;
		uint16_t* B_i = pk_CPA.B + LWE_L * i;
		for(int k = 0; k < LWE_N; ++k){
			uint16_t* sk_k = sk_CPA + LWE_L * k;
			uint16_t A_ik = A_i[k];
			for(int j = 0; j < LWE_L; ++j){
				B_i[j] -= A_ik * sk_k[j];
			}
		}
	}
	}
	elapsed1 = clock() - elapsed1;
    printf("    Keygen Time: %f ms\n",elapsed1 * 1000./CLOCKS_PER_SEC/iter);	
}

// generate LWE_M dimensional coefficient vector r
// store indices of nonzero components of r in HR dimensional vector r_idx
// for i in [0, neg_start), the vector r has the value 1 at index r_idx[i] 
// for i in [neg_start, HR), the vector r has the value 1 at index r_idx[i] 

size_t gen_r_idx(uint16_t* r_idx){
	size_t neg_start = 0;
	for(int i = 0; i < HR; ++i){
		r_idx[i] = rand() % LWE_M;
		neg_start += (rand() & 0x01);
	}
	return neg_start;
}

void Enc_CPA() {
	uint16_t r_idx[HR];
	size_t neg_start;

	start = clock();
	// Initialize ct.a as zero vector, ct.b as q/t * m.
	for(int i = 0; i < LWE_N; ++i) { ctx_CPA.a[i] = 0; }
	for(int i = 0; i < LWE_L; ++i) { ctx_CPA.b[i] = ptx_CPA[i] << _16_LOG_T;}

	neg_start = gen_r_idx(r_idx);
	
	// Compute A^T * r and B^T * r, and then add to ct.a and ct.b, respectively.
	// use pointers for optimization of matrix multiplication
	for (size_t i = 0; i < neg_start; ++i) {
		uint16_t* pk_A_ri = pk_CPA.A + LWE_N * r_idx[i];
		uint16_t* pk_B_ri = pk_CPA.B + LWE_L * r_idx[i];
		for (int j = 0; j < LWE_N; ++j) { ctx_CPA.a[j] += pk_A_ri[j]; }
		for (int j = 0; j < LWE_L; ++j) { ctx_CPA.b[j] += pk_B_ri[j]; }
	}

	for (size_t i = neg_start; i < HR; ++i) {
		uint16_t* pk_A_ri = pk_CPA.A + LWE_N * r_idx[i];
		uint16_t* pk_B_ri = pk_CPA.B + LWE_L * r_idx[i];
		for (int j = 0; j < LWE_N; ++j) { ctx_CPA.a[j] -= pk_A_ri[j]; }
		for (int j = 0; j < LWE_L; ++j) { ctx_CPA.b[j] -= pk_B_ri[j]; }
	}

	// Send a and b from mod q to mod p
	for(int i = 0; i < LWE_N; ++i) { ctx_CPA.a[i] = (ctx_CPA.a[i] + RD_ADD) & RD_AND; }
	for(int i = 0; i < LWE_L; ++i) { ctx_CPA.b[i] = (ctx_CPA.b[i] + RD_ADD) & RD_AND; }

	finish = clock();
	elapsed1 += (finish-start);
}

void Dec_CPA(){
	start = clock();
	for(int i = 0; i < LWE_L; ++i) { ptx_CPA[i] = ctx_CPA.b[i]; } // initialize msg as ct.b	

	for(int i = 0; i < LWE_N; ++i){
		uint16_t* sk_i = sk_CPA + LWE_L * i;
		uint16_t ctx_ai = ctx_CPA.a[i];
		for(int j = 0; j < LWE_L; ++j){
			ptx_CPA[j] += ctx_ai * sk_i[j];
		}
	}

	for(int i = 0; i < LWE_L; ++i){ 
		ptx_CPA[i] += DEC_ADD; 
		ptx_CPA[i] >>= _16_LOG_T; 
	}

	finish = clock();
	elapsed2 += (finish-start);
}

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////   CCA SCHEME   ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

uint64_t ptx_CCA[4];

struct { // n-dimensional LWE ctxtext modulo q
    uint16_t b[LWE_L];    // integer modulo q
    uint16_t a[LWE_N]; // array of n integers modulo q
    uint64_t c1[4];
    uint64_t c3[6];
} ctx_CCA;

uint64_t backup_CCA[4];

int mod (int a, int b){
	if(b < 0)
		return mod(-a, -b);
	int ret = a % b;
	if(ret < 0)
		ret += b;
	return ret;
}

const uint64_t r[25]={ 0,36,3,41,18, 1,44,10,45,2, 62,6,43,15,61, 28,55,25,21,56, 27,20,39,8,14 };

//Round constant in keccak//
const uint64_t RC[24]={ 0x0000000000000001, 0x0000000000008082, 0x800000000000808A, 0x8000000080008000,
		       			0x000000000000808B, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
		       			0x000000000000008A, 0x0000000000000088, 0x0000000080008009, 0x000000008000000A,
		       			0x000000008000808B, 0x800000000000008B, 0x8000000000008089, 0x8000000000008003,
		       			0x8000000000008002, 0x8000000000000080, 0x000000000000800A, 0x800000008000000A,
		       			0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008};

uint64_t* sha3_round(uint64_t* A, uint64_t RC){
	uint8_t x, y;
	uint64_t C[5];
	uint64_t D[5];
	uint64_t B[25];

	//Theta step
	for(x = 0; x < 5; x++){C[x] = A[5 * x] ^ A[5 * x + 1] ^ A[5 * x + 2] ^ A[5 * x + 3] ^ A[5 * x + 4];}
	for(x = 0; x < 5; x++){D[x] = C[(x + 4) % 5] ^ ((C[(x + 1) % 5] << 1) | (C[(x + 1) % 5] >> 63));}
	for(x = 0; x < 5; x++){ for(y = 0; y < 5; y++){ A[5 * x + y] = A[5 * x + y] ^ D[x]; }}

	//Rho and phi step
	for(x = 0; x < 5; x++){ for(y = 0; y < 5; y++){
		B[5 * y + mod((2 * x + 3 * y),5)] = ((A[5 * x + y] << r[5 * x + y]) | (A[5 * x + y] >> (64 - r[5 * x + y])));
	}}

	//Xi step
	for(x = 0; x < 5; x++){ for(y = 0; y < 5; y++ ){
		A[5 * x + y] = B[5 * x + y] ^ ((~B[5 * mod((x+1),5) + y]) & B[5 * mod((x+2),5) + y]);
	}}

	//XOR step
	A[0] = A[0] ^ RC;
	return A;
}

//keccak fuction with sha3 round function: input and output is 5 by 5 uint64 matrix
uint64_t *keccak_f(uint64_t *A){
	for(int32_t i = 0; i < 24; i++){ A = sha3_round(A, RC[i]); }
	return A;
}

void Enc_CCA() {
	start = clock();
	// random assign
	uint64_t delta[6] = {rand(), rand(), rand(), rand(), rand(), rand()};
  	uint64_t *hash = (uint64_t *)calloc(25, sizeof(uint64_t)); //Set seed to be 5 by 5 matrix of uint64_t

	memset((unsigned char*) hash, 0, 200);
	memcpy((unsigned char*) hash, (unsigned char*) delta, 48);
    
	hash = keccak_f(hash);
	ctx_CCA.c1[0] = ptx_CCA[0] ^ hash[0];
	ctx_CCA.c1[1] = ptx_CCA[1] ^ hash[1];
	ctx_CCA.c1[2] = ptx_CCA[2] ^ hash[2];
	ctx_CCA.c1[3] = ptx_CCA[3] ^ hash[3];

	ctx_CCA.c3[0] = hash[4];
	ctx_CCA.c3[1] = hash[5];
	ctx_CCA.c3[2] = hash[6];
	ctx_CCA.c3[3] = hash[7];
	ctx_CCA.c3[4] = hash[8];
	ctx_CCA.c3[5] = hash[9];

	memset((unsigned char*) hash, 0, 200);
	memcpy((unsigned char*) hash, (unsigned char*) delta, 48);
	memcpy(((unsigned char*) hash) + 48, (unsigned char*) ctx_CCA.c1, 32);

	hash = keccak_f(hash);
		
	// Initialize ct.a as zero vector, ct.b as q/t * m.
	for(int j = 0; j < LWE_N; ++j) { ctx_CCA.a[j] = 0; }
	for(int i = 0; i < 6; ++i){ for(int j = 0; j < 64; ++j){ ctx_CCA.b[64 * i + j] = ((uint16_t) (delta[i] >> j)) << 15; }}
	
	// decomp_delta[64 * i + j] = (delta[i] >> j) & 0x01;

	uint16_t r_idx[HR];
	size_t neg_start = 0;
	
	for(int i = 0; i < 21; ++i){
		int j = 6 * i;
		r_idx[j] = (uint16_t) hash[i] & 0x03ff;
		r_idx[j+1] = (uint16_t) (hash[i] >> 10) & 0x03ff;
		r_idx[j+2] = (uint16_t) (hash[i] >> 20) & 0x03ff;
		r_idx[j+3] = (uint16_t) (hash[i] >> 30) & 0x03ff;
		r_idx[j+4] = (uint16_t) (hash[i] >> 40) & 0x03ff;
		r_idx[j+5] = (uint16_t) (hash[i] >> 50) & 0x03ff;
	}
	r_idx[126] = (uint16_t) hash[21] & 0x03ff;
	r_idx[127] = (uint16_t) (hash[21] >> 10) & 0x03ff;

	for(int i = 0; i < 64; ++i){
		neg_start += (size_t) ((hash[22] >> i) & 0x01);
		neg_start += (size_t) ((hash[23] >> i) & 0x01);
	}

	//cout << "Hash end" << endl;
	// Compute A^T * r and B^T * r, and then add to ct.a and ct.b, respectively.
	for (size_t i = 0; i < neg_start; ++i) {
		uint16_t idx = r_idx[i];
		for (int j = 0; j < LWE_N; ++j) { ctx_CCA.a[j] += pk_CPA.A[LWE_N * idx + j]; }
		for (int j = 0; j < LWE_L; ++j) { ctx_CCA.b[j] += pk_CPA.B[LWE_L * idx + j]; }
	}
	for (size_t i = neg_start; i < HR; ++i) {
		uint16_t idx = r_idx[i];
		for (int j = 0; j < LWE_N; ++j) { ctx_CCA.a[j] -= pk_CPA.A[LWE_N * idx + j]; }
		for (int j = 0; j < LWE_L; ++j) { ctx_CCA.b[j] -= pk_CPA.B[LWE_L * idx + j]; }
	}
	// Send a and b from mod q to mod p
	for(int i = 0; i < LWE_N; ++i) { ctx_CCA.a[i] = (ctx_CCA.a[i] + RD_ADD) & RD_AND; }
	for(int i = 0; i < LWE_L; ++i) { ctx_CCA.b[i] = (ctx_CCA.b[i] + RD_ADD) & RD_AND; }
	finish = clock();
	elapsed1 += (finish-start);
	free(hash);
}
	
int Dec_CCA(){
	uint16_t decomp_delta[LWE_L];	
	int res = 0;

	start = clock();

	for(int i = 0; i < LWE_L; ++i) { decomp_delta[i] = ctx_CCA.b[i]; } // initialize msg as ct.b	
	for(int i = 0; i < LWE_N; ++i){
		for(int j = 0; j < LWE_L; ++j){
			decomp_delta[j] += ctx_CCA.a[i] * sk_CPA[LWE_L * i + j];
		}
	}

	for(int i = 0; i < LWE_L; ++i){ 
		decomp_delta[i] +=  DEC_ADD; 
		decomp_delta[i] >>= _16_LOG_T; 
	}


	uint64_t delta[6] = {};
	for(int i = 0; i < 6; ++i){ for(int j = 0; j < 64; ++j){
		uint64_t a = ((uint64_t) decomp_delta[64 * i + j]) << j;
		delta[i] ^= a;
	}}
	// decomp_delta[64 * i + j] = (delta[i] >> j) & 0x01;

  	uint64_t *hash = (uint64_t *)calloc(25, sizeof(uint64_t)); //Set seed to be 5 by 5 matrix of uint64_t

	memset((unsigned char*) hash, 0, 200);
	memcpy((unsigned char*) hash, (unsigned char*) delta, 48);

	hash = keccak_f(hash);

	for(int i = 0; i < 4; ++i){ ptx_CCA[i] = ctx_CCA.c1[i] ^ hash[i]; }

	if((hash[4] != ctx_CCA.c3[0]) || (hash[5] != ctx_CCA.c3[1]) || (hash[6] != ctx_CCA.c3[2])|| 
		(hash[7] != ctx_CCA.c3[3])|| (hash[8] != ctx_CCA.c3[4])|| (hash[9] != ctx_CCA.c3[5])){ res =  1; }

	memset((unsigned char*) hash, 0, 200);
	memcpy((unsigned char*) hash, (unsigned char*) delta, 48);
	memcpy(((unsigned char*) hash) + 48, (unsigned char*) ctx_CCA.c1, 32);

	hash = keccak_f(hash);

	uint16_t c2h_b[LWE_L];
	for(int i = 0; i < LWE_L; ++i){ c2h_b[i] = decomp_delta[i] << 15; }
	uint16_t c2h_a[LWE_N] = {};

	uint16_t r_idx[128];
	size_t neg_start = 0;
	
	for(int i = 0; i < 21; ++i){
		int j = 6 * i;
		r_idx[j] = (uint16_t) hash[i] & 0x03ff;
		r_idx[j+1] = (uint16_t) (hash[i] >> 10) & 0x03ff;
		r_idx[j+2] = (uint16_t) (hash[i] >> 20) & 0x03ff;
		r_idx[j+3] = (uint16_t) (hash[i] >> 30) & 0x03ff;
		r_idx[j+4] = (uint16_t) (hash[i] >> 40) & 0x03ff;
		r_idx[j+5] = (uint16_t) (hash[i] >> 50) & 0x03ff;
	}
	r_idx[126] = (uint16_t) hash[21] & 0x03ff;
	r_idx[127] = (uint16_t) (hash[21] >> 10) & 0x03ff;

	for(int i = 0; i < 64; ++i){
		neg_start += (size_t) ((hash[22] >> i) & 0x01);
		neg_start += (size_t) ((hash[23] >> i) & 0x01);
	}
	
	//cout << "Hash end" << endl;
	// Compute A^T * r and B^T * r, and then add to ct.a and ct.b, respectively.
	for (int i = 0; i < (int)neg_start; ++i) {
		uint16_t idx = r_idx[i];
		for (int j = 0; j < LWE_N; ++j) { c2h_a[j] += pk_CPA.A[LWE_N * idx + j]; }
		for (int j = 0; j < LWE_L; ++j) { c2h_b[j] += pk_CPA.B[LWE_L * idx + j]; }
	}
	for (int i = neg_start; i < 128; ++i) {
		uint16_t idx = r_idx[i];
		for (int j = 0; j < LWE_N; ++j) { c2h_a[j] -= pk_CPA.A[LWE_N * idx + j]; }
		for (int j = 0; j < LWE_L; ++j) { c2h_b[j] -= pk_CPA.B[LWE_L * idx + j]; }
	}

	// Send a and b from mod q to mod p
	for(int i = 0; i < LWE_N; ++i) { c2h_a[i] = (c2h_a[i] + RD_ADD) & RD_AND; }
	for(int i = 0; i < LWE_L; ++i) { c2h_b[i] = (c2h_b[i] + RD_ADD) & RD_AND; }

	for(int i = 0; i < LWE_N; ++i){ if(c2h_a[i] != ctx_CCA.a[i]) res =  2; }
	for(int i = 0; i < LWE_L; ++i){ if(c2h_b[i] != ctx_CCA.b[i]) res =  2; }

	if(res != 0) { for(int i = 0; i < 4; ++i){ ptx_CCA[i] = 0; }}
	finish = clock();
	elapsed2 += (finish-start);

	free(hash);
	return res;
}

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////   RING LIZARD SCHEME   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

struct {
	uint16_t idx[HR];
	size_t neg_start;
} sk_ring;

struct {
	uint16_t a[LWE_N];
	uint16_t b[LWE_N];
} pk_ring;

uint16_t ptx_ring[LWE_N];

uint16_t backup_ring[LWE_N];

struct { // n-dimensional LWE ctxtext modulo q
	uint16_t a[LWE_N];		// array of n integers modulo q
	uint16_t b[LWE_N];		// integer modulo q
} ctx_ring;


void gen_sk_ring(){
	sk_ring.neg_start = 0;
	for(int i = 0; i < HR; ++i){
		sk_ring.idx[i] = rand() % LWE_N;
		sk_ring.neg_start += rand() & 0x01;
	}
}

void gen_a_ring(){
	for(int i = 0; i < LWE_N; ++i){
		pk_ring.a[i] = rand() << _16_LOG_Q;
	}
}

void gen_e_ring(){
	for(int i = 0; i < LWE_N; ++i){
		pk_ring.b[i] = SAMPLE_DG() << _16_LOG_Q;
	}
}

void Keygen_ring(){
	elapsed1 = clock();
	for(int l = 0; l < iter; ++l){
	gen_sk_ring();
	gen_a_ring();
	gen_e_ring();

	for (size_t i = 0; i < sk_ring.neg_start; ++i) {
		uint16_t deg = sk_ring.idx[i];
		for (int j = 0; j < LWE_N - deg; ++j) { pk_ring.b[deg + j] -= pk_ring.a[j]; }
		for (int j = LWE_N - deg; j < LWE_N; ++j){ pk_ring.b[deg + j - LWE_N] += pk_ring.a[j]; }
	}

	for (size_t i = sk_ring.neg_start; i < HR; ++i) {
		uint16_t deg = sk_ring.idx[i];
		for (int j = 0; j < LWE_N - deg; ++j) { pk_ring.b[deg + j] += pk_ring.a[j]; }
		for (int j = LWE_N - deg; j < LWE_N; ++j){ pk_ring.b[deg + j - LWE_N] -= pk_ring.a[j]; }
	}
	}

	elapsed1 = clock() - elapsed1;
	printf("    Keygen Time: %f ms\n",elapsed1 * 1000./CLOCKS_PER_SEC/iter);	
}


void Enc_ring() {
	uint16_t r_idx[HR];
	size_t neg_start;

	start = clock();
	// Initialize ct.a as zero vector, ct.b as q/t * m.
	for(int i = 0; i < LWE_N; ++i) { ctx_ring.a[i] = 0; }
	for(int i = 0; i < LWE_N; ++i) { ctx_ring.b[i] = ptx_ring[i] << _16_LOG_T;}

	neg_start = gen_r_idx(r_idx);
	
	// Compute a * r and b * r, and then add to ct.a and ct.b, respectively.

	for (size_t i = 0; i < neg_start; ++i) {
		uint16_t deg = r_idx[i];
		for (int j = 0; j < LWE_N - deg; ++j) { ctx_ring.a[deg + j] += pk_ring.a[j]; }
		for (int j = LWE_N - deg; j < LWE_N; ++j){ ctx_ring.a[deg + j - LWE_N] -= pk_ring.a[j]; }

		for (int j = 0; j < LWE_N - deg; ++j) { ctx_ring.b[deg + j] += pk_ring.b[j]; }
		for (int j = LWE_N - deg; j < LWE_N; ++j){ ctx_ring.b[deg + j - LWE_N] -= pk_ring.b[j]; }
	}

	for (size_t i = neg_start; i < HR; ++i) {
		uint16_t deg = r_idx[i];
		for (int j = 0; j < LWE_N - deg; ++j) { ctx_ring.a[deg + j] -= pk_ring.a[j]; }
		for (int j = LWE_N - deg; j < LWE_N; ++j){ ctx_ring.a[deg + j - LWE_N] += pk_ring.a[j]; }

		for (int j = 0; j < LWE_N - deg; ++j) { ctx_ring.b[deg + j] -= pk_ring.b[j]; }
		for (int j = LWE_N - deg; j < LWE_N; ++j){ ctx_ring.b[deg + j - LWE_N] += pk_ring.b[j]; }
	}

	// Send a and b from mod q to mod p
	for(int i = 0; i < LWE_N; ++i) { ctx_ring.a[i] = (ctx_ring.a[i] + RD_ADD) & RD_AND; }
	for(int i = 0; i < LWE_N; ++i) { ctx_ring.b[i] = (ctx_ring.b[i] + RD_ADD) & RD_AND; }

	finish = clock();
	elapsed1 += (finish-start);
}

void Dec_ring(){
	start = clock();
	for(int i = 0; i < LWE_N; ++i){ ptx_ring[i] = ctx_ring.b[i]; }

	for (size_t i = 0; i < sk_ring.neg_start; ++i) {
		uint16_t deg = sk_ring.idx[i];
		for (int j = 0; j < LWE_N - deg; ++j) { ptx_ring[deg + j] += ctx_ring.a[j]; }
		for (int j = LWE_N - deg; j < LWE_N; ++j){ ptx_ring[deg + j - LWE_N] -= ctx_ring.a[j]; }
	}

	for (size_t i = sk_ring.neg_start; i < HR; ++i) {
		uint16_t deg = sk_ring.idx[i];
		for (int j = 0; j < LWE_N - deg; ++j) { ptx_ring[deg + j] -= ctx_ring.a[j]; }
		for (int j = LWE_N - deg; j < LWE_N; ++j){ ptx_ring[deg + j - LWE_N] += ctx_ring.a[j]; }
	}

	for(int i = 0; i < LWE_N; ++i){ 
		ptx_ring[i] += DEC_ADD; 
		ptx_ring[i] >>= _16_LOG_T; 
	}

	finish = clock();
	elapsed2 += (finish-start);
}

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////   TEST FUNCTIONS   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void EncDecTest_CPA(){
    int i;

    elapsed1 = 0;
    elapsed2 = 0;

    for(int l = 0; l < iter; ++l){
    // Generate random messages
    for(i = 0; i < LWE_L; ++i) {
        ptx_CPA[i] = rand() % T;
		backup_CPA[i] = ptx_CPA[i];
    }

    for( i = 0; i < testnum; i++ ) Enc_CPA();
    
    for( i = 0; i < testnum; i++ ) Dec_CPA();

    // Correctness check
  
    for(i = 0; i < LWE_L; ++i) {
        if(ptx_CPA[i] != backup_CPA[i]) {
	    printf("    Correctness Error\n");
            break;
        }
    }
    if(i < LWE_L) break;
	}
	printf("    Enc Time: %f ms\n",elapsed1 * 1000./CLOCKS_PER_SEC/testnum/iter);
    printf("    Dec Time: %f ms\n",elapsed2 * 1000./CLOCKS_PER_SEC/testnum/iter);
}

void EncDecTest_CCA(){
	int i, res;
    
    // Generate random messages
    for(i = 0; i < 4; ++i) {
        ptx_CCA[i] = rand();
		backup_CCA[i] = ptx_CCA[i];
    }

    elapsed1 = 0;
    elapsed2 = 0;

    for(int l = 0; l < iter; ++l){

    // Generate random messages
    for( i = 0; i < testnum; i++ ) Enc_CCA();
    
    for( i = 0; i < testnum; i++ ) res = Dec_CCA();

	// Correctness check

	if(res == 1 ){
		printf("    Decryption Validity Error Type 1 : c3 components\n");
		break;
	}
	if(res == 2 ){
		printf("    Decryption Validity Error Type 2 : c2 components\n");
		break;
	}

	for(i = 0; i < 4; ++i) {
		if(ptx_CCA[i] != backup_CCA[i]) {
			printf("    Correctness Error\n");
			break;
		}
    }
    if(i < 4) break;

	}
	printf("    Enc Time: %f ms\n", elapsed1 * 1000./CLOCKS_PER_SEC/testnum/iter);
    printf("    Dec Time: %f ms\n", elapsed2 * 1000./CLOCKS_PER_SEC/testnum/iter);
}

void EncDecTest_RING(){
	int i;

	elapsed1 = 0;
	elapsed2 = 0;

	for(int l = 0; l < iter; ++l){
	// Generate random messages
	
	for(i = 0; i < LWE_N; ++i) {
        ptx_ring[i] = rand() % T;
		backup_ring[i] = ptx_ring[i];
    }

	for( i = 0; i < testnum; i++ ) Enc_ring();
	for( i = 0; i < testnum; i++ ) Dec_ring();

    // Correctness check
  
    for(i = 0; i < LWE_N; ++i) {
        if(ptx_ring[i] != backup_ring[i]) {
	    printf("Correctness Error, %d\n",i);
            break;
        }
    }
    if(i < LWE_N) break;
	}

	printf("    Enc Time: %f ms\n",elapsed1 * 1000./CLOCKS_PER_SEC/testnum/iter);
    printf("    Dec Time: %f ms\n",elapsed2 * 1000./CLOCKS_PER_SEC/testnum/iter);
}

#ifdef PARAMS_Classical
#define KEYGEN Keygen_CPA
#define ENCDECTEST EncDecTest_CPA
#endif

#ifdef PARAMS_Recommended
#define KEYGEN Keygen_CPA
#define ENCDECTEST EncDecTest_CPA
#endif

#ifdef PARAMS_Classical_Plaintext_32bit
#define KEYGEN Keygen_CPA
#define ENCDECTEST EncDecTest_CPA
#endif

#ifdef PARAMS_CCA
#define KEYGEN Keygen_CPA
#define ENCDECTEST EncDecTest_CCA
#endif

#ifdef PARAMS_RING
#define KEYGEN Keygen_ring
#define ENCDECTEST EncDecTest_RING
#endif

////////////////////////////////////////////////////////////////////////////////////////////

int main(){
	srand((unsigned)time(NULL)+(unsigned)getpid());
	printf("\n  //////////////////////////////////////////////////////////////////\n\n");
	printf("\t\t\t"PARAMNAME" Parameter\n\n");
	printf("    LWE dimension: %d, \t\tLWR dimension: %d\n", LWE_N, LWE_M);
	printf("    Plaintext dimension: %d, \t\tPlaintext Modulus: %d bits\t\n", LWE_N, LOG_T);
	printf("    Public Key modulus: %d bits, \tCiphertext modulus: %d bits\t\n\n", LOG_Q, LOG_P);
	printf("  //////////////////////////////////////////////////////////////////\n\n");
	printf("\t\t\tPerformance Test\n\n");

	// Key Generation
	KEYGEN();

	// Encryption and Decryption Test
	ENCDECTEST();	
	printf("\n  //////////////////////////////////////////////////////////////////\n\n");

}
