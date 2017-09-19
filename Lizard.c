#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <immintrin.h>
#include <stdint.h>
#include <inttypes.h>
#include "fips202.h"


// compile and run command with Makefile:
// $ make new


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// CYCLE COUNTING MODULE ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

long long rdtsc(void)
{
  unsigned long long result;
  __asm__ volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
    : "=a" (result) ::  "%rdx");
  return result;
}


//////////////////////////////////  PARAMETER CHOICES  /////////////////////////////////////
//#define PARAMS_Classical								// classical parameter for Lizard
#define PARAMS_Recommended								// recommended parameter for Lizard
//#define PARAMS_Homadd									// recommended parameter for Lizard which supports 100 number of homomorphic additions
//#define PARAMS_Classical_Plaintext_32bit				// classical parameter for Lizard with 32-bit plaintext
//#define PARAMS_CCA									// recommended parameter for CCALizard


//////////////////////////////////  PARAMETER SETS  ////////////////////////////////////////

#define iter 100		// iteration number for keygen & EncDec test
#define testnum 1000	// repeatetion number of Enc Dec procedure in a single iteration

#ifdef PARAMS_Classical
#define PARAMNAME "Classical"
#define LWE_N 480		// LWE dim
#define LWE_M 724		// LWR dim (Number of LWE samples in pk)
#define LWE_L 256
#define LOG_Q 11
#define _16_LOG_Q 5
#define LOG_P 9
#define RD_ADD 0x40 	// 2^(15 - LOG_P)
#define RD_AND 0xff80
#define LOG_T 1
#define _16_LOG_T 15
#define T 2
#define DEC_ADD 0x4000	// 2^(15 - LOG_T)
#define NOISE_D1		// standard deviation for discrete gaussian distribution
#define HR 128			// Hamming weight of coefficient vector r
#endif

#ifdef PARAMS_Recommended
#define PARAMNAME "Recommended"
#define LWE_N 536		// LWE dim
#define LWE_M 1024		// LWR dim (Number of LWE samples in pk)
#define LWE_L 256
#define LOG_Q 11
#define _16_LOG_Q 5
#define LOG_P 9
#define RD_ADD 0x40 	// 2^(15 - LOG_P)
#define RD_AND 0xff80
#define LOG_T 1
#define _16_LOG_T 15
#define T 2
#define DEC_ADD 0x4000	// 2^(15 - LOG_T)
#define NOISE_D2		// standard deviation for discrete gaussian distribution
#define HR 134			// Hamming weight of coefficient vector r
#endif

#ifdef PARAMS_Homadd
#define PARAMNAME "Homadd"
#define LWE_N 816		// LWE dim
#define LWE_M 1024		// LWR dim (Number of LWE samples in pk)
#define LWE_L 256
#define LOG_Q 16
#define _16_LOG_Q 0
#define LOG_P 14
#define RD_ADD 0x02 	// 2^(15 - LOG_P)
#define RD_AND 0xfffc
#define LOG_T 1
#define _16_LOG_T 15
#define T 2
#define DEC_ADD 0x4000	// 2^(15 - LOG_T)
#define NOISE_D2		// standard deviation for discrete gaussian distribution
#define HR 136			// Hamming weight of coefficient vector r
#endif


#ifdef PARAMS_Classical_Plaintext_32bit
#define PARAMNAME "Classical_32bit"
#define LWE_N 480		// LWE dim
#define LWE_M 724		// LWR dim (Number of LWE samples in pk)
#define LWE_L 32
#define LOG_Q 11
#define _16_LOG_Q 5
#define LOG_P 9
#define RD_ADD 0x40 	// 2^(15 - LOG_P)
#define RD_AND 0xff80
#define LOG_T 1
#define _16_LOG_T 15
#define T 2
#define DEC_ADD 0x4000
#define NOISE_D1		// standard deviation for discrete gaussian distribution
#define HR 128			// Hamming weight of coefficient vector r
#endif

#ifdef PARAMS_CCA
#define PARAMNAME "CCA"
#define LWE_N 536		// LWE dim
#define LWE_M 1024		// LWR dim (Number of LWE samples in pk)
#define LWE_L 256
#define LOG_Q 11
#define _16_LOG_Q 5
#define LOG_P 9
#define RD_ADD 0x40 	// 2^(15 - LOG_P)
#define RD_AND 0xff80
#define HR 134
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

//A discrete error distribution close to the discrete Gaussian distribution with sigma = 6.759075908 w.r.t. Renyi div = 1.00327, Renyi deg = 25
#ifdef NOISE_D1
#define SAMPLE_DG Sample_D1
const uint16_t CDF_TABLE[9] = {75, 216, 332, 414, 465, 492, 504, 509, 511}; // out of [0, 511]
const size_t RANDBITS = 10;
const size_t TABLE_LENGTH = 9;

uint16_t Sample_D1(){
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

//A discrete error distribution close to the discrete Gaussian distribution with sigma = 6.481012658 w.r.t. Renyi div = 1.005, Renyi deg = 25
#ifdef NOISE_D2
#define SAMPLE_DG Sample_D2
const uint16_t CDF_TABLE[9] = {78, 226, 344, 425, 473, 495, 506, 510, 511}; // out of [0, 511]
const size_t RANDBITS = 10;
const size_t TABLE_LENGTH = 9;

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

//A discrete error distribution close to the discrete Gaussian distribution with sigma = 3.357377049 w.r.t. Renyi div = 1.005624513, Renyi deg = 25
#ifdef NOISE_D3
#define SAMPLE_DG Sample_D3
const uint16_t CDF_TABLE[5] = {151, 382, 482, 507, 511}; // out of [0, 511]
const size_t RANDBITS = 10;
const size_t TABLE_LENGTH = 5;

uint16_t Sample_D3(){
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

//A discrete error distribution close to the discrete Gaussian distribution with sigma = 7.801904762 w.r.t. Renyi div = 1.001424977, Renyi deg = 25
#ifdef NOISE_D4
#define SAMPLE_DG Sample_D4
const uint16_t CDF_TABLE[11] = {131, 380, 594, 759, 874, 946, 987, 1008, 1018, 1022, 1023}; // out of [0, 1023]
const size_t RANDBITS = 11;
const size_t TABLE_LENGTH = 11;

uint16_t Sample_D4(){
	uint16_t rnd = rand() & 0x03ff;
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
uint64_t start_cycle, finish_cycle, cycles1, cycles2;

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// CPA SCHEME //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

// Every variables are defined over uint16_t
// For smaller modulus operations, fill the zeros in LSBs of uint16_t
// e.g. when q=2^10, pad 6 bits in LSBs during pk generation

uint16_t sk_CPA[LWE_N * LWE_L];

uint64_t s[4];

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
    cycles1 = rdtsc();
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
	cycles1 = rdtsc() - cycles1;
    printf("    Keygen Time: %f ms\n",elapsed1 * 1000./CLOCKS_PER_SEC/iter);	
    printf("    Keygen Cycles:  %llu \n",  cycles1/iter);	
}


// generate LWE_M dimensional coefficient vector r
// store indices of nonzero components of r in HR dimensional vector r_idx
// for i in [0, neg_start), the vector r has the value 1 at index r_idx[i] 
// for i in [neg_start, HR), the vector r has the value 1 at index r_idx[i] 

size_t gen_r_idx(uint16_t* r_idx){
	size_t neg_start = 0;
	uint64_t tmp;

	for(int i = 0; i < HR/2; ++i){
		tmp = rand();
		neg_start += (tmp & 0x01);
		r_idx[2 * i] = ((tmp >> 1) & 0x03ff) % LWE_M;
		neg_start += ((tmp >> 0x7ff) & 0x01);
		r_idx[2 * i + 1] = ((tmp >> 0xfff) & 0x03ff) % LWE_M;
	}
	return neg_start;
}

void Enc_CPA() {
	start = clock();
	start_cycle = rdtsc();

	uint16_t r_idx[HR];
	size_t neg_start;
	
	// Initialize ct.a as zero vector, ct.b as q/t * m.
	for(int i = 0; i < LWE_N; ++i) { ctx_CPA.a[i] = 0; }
	for(int i = 0; i < LWE_L; ++i) { ctx_CPA.b[i] = ptx_CPA[i] << _16_LOG_T;}

	neg_start = gen_r_idx(r_idx);

	// Compute A^T * r and B^T * r, and then add to ct.a and ct.b, respectively.
	// use pointers for optimization of matrix multiplication

	for (size_t i = 0; i < HR; ++i) {
        uint16_t  s = (i < neg_start)?1:0;
        uint16_t* pk_A_ri = pk_CPA.A + LWE_N * r_idx[i];
        uint16_t* pk_B_ri = pk_CPA.B + LWE_L * r_idx[i];
        for (int j = 0; j < LWE_N; ++j) { 
                ctx_CPA.a[j] += s * pk_A_ri[j];
                ctx_CPA.a[j] -= (1 - s) * pk_A_ri[j];
        }
        for (int j = 0; j < LWE_L; ++j) { 
                ctx_CPA.b[j] += s * pk_B_ri[j];
                ctx_CPA.b[j] -= (1 - s) * pk_B_ri[j];
        }
	}

	// Send a and b from mod q to mod p
	for(int i = 0; i < LWE_N; ++i) { ctx_CPA.a[i] = (ctx_CPA.a[i] + RD_ADD) & RD_AND; }
	for(int i = 0; i < LWE_L; ++i) { ctx_CPA.b[i] = (ctx_CPA.b[i] + RD_ADD) & RD_AND; }

	finish = clock();
	finish_cycle = rdtsc();
	elapsed1 += (finish-start);
	cycles1 += (finish_cycle - start_cycle);
}

void Dec_CPA(){
	start = clock();
	start_cycle = rdtsc();
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
	finish_cycle = rdtsc();
	elapsed2 += (finish-start);
	cycles2 += (finish_cycle - start_cycle);
}

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////   CCA SCHEME   ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

uint64_t ptx_CCA[4];

struct { // n-dimensional LWE ctxtext modulo q
    uint16_t b[LWE_L];    // integer modulo q
    uint16_t a[LWE_N];	  // array of n integers modulo q
    uint64_t c1[4];
    uint64_t c3[4];
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

void Enc_CCA() {
	start = clock();
	start_cycle = rdtsc();
	// random assign
	uint64_t delta[4] = {rand(), rand(), rand(), rand()};
  	uint64_t *hash = (uint64_t *)calloc(32, sizeof(uint64_t)); 

	memset((unsigned char*) hash, 0, 256);
	memcpy((unsigned char*) hash, (unsigned char*) delta, 32);
    shake128((unsigned char*)hash, 256, (unsigned char*)hash, 32);

	ctx_CCA.c1[0] = ptx_CCA[0] ^ hash[0];
	ctx_CCA.c1[1] = ptx_CCA[1] ^ hash[1];
	ctx_CCA.c1[2] = ptx_CCA[2] ^ hash[2];
	ctx_CCA.c1[3] = ptx_CCA[3] ^ hash[3];

	ctx_CCA.c3[0] = hash[4];
	ctx_CCA.c3[1] = hash[5];
	ctx_CCA.c3[2] = hash[6];
	ctx_CCA.c3[3] = hash[7];

		
	// Initialize ct.a as zero vector, ct.b as q/t * m.
	for(int j = 0; j < LWE_N; ++j) { ctx_CCA.a[j] = 0; }
	for(int i = 0; i < 4; ++i){ for(int j = 0; j < 64; ++j){ ctx_CCA.b[64 * i + j] = ((uint16_t) (delta[i] >> j)) << 15; }}
	
	// decomp_delta[64 * i + j] = (delta[i] >> j) & 0x01;

	uint16_t r_idx[HR];
	size_t neg_start = 0;
	
	for(int i = 8; i < 30; ++i){
		int j = 6 * (i - 8);
		r_idx[j] = (uint16_t) hash[i] & 0x03ff;
		r_idx[j+1] = (uint16_t) (hash[i] >> 10) & 0x03ff;
		r_idx[j+2] = (uint16_t) (hash[i] >> 20) & 0x03ff;
		r_idx[j+3] = (uint16_t) (hash[i] >> 30) & 0x03ff;
		r_idx[j+4] = (uint16_t) (hash[i] >> 40) & 0x03ff;
		r_idx[j+5] = (uint16_t) (hash[i] >> 50) & 0x03ff;
		neg_start += (size_t) ((hash[i] >> 60) & 0x01);
		neg_start += (size_t) ((hash[i] >> 61) & 0x01);
		neg_start += (size_t) ((hash[i] >> 62) & 0x01);
		neg_start += (size_t) ((hash[i] >> 63) & 0x01);
	}
	r_idx[132] = (uint16_t) hash[30] & 0x03ff;
	r_idx[133] = (uint16_t) (hash[30] >> 10) & 0x03ff;


	for(int i = 0; i < 46; ++i){
		neg_start += (size_t) ((hash[31] >> i) & 0x01);
	}

	//cout << "Hash end" << endl;
	// Compute A^T * r and B^T * r, and then add to ct.a and ct.b, respectively.
	for (size_t i = 0; i < HR; ++i) {
        uint16_t  s = (i < neg_start)?1:0;
        uint16_t* pk_A_ri = pk_CPA.A + LWE_N * r_idx[i];
        uint16_t* pk_B_ri = pk_CPA.B + LWE_L * r_idx[i];
        for (int j = 0; j < LWE_N; ++j) { 
                ctx_CCA.a[j] += s * pk_A_ri[j];
                ctx_CCA.a[j] -= (1 - s) * pk_A_ri[j];
        }
        for (int j = 0; j < LWE_L; ++j) { 
                ctx_CCA.b[j] += s * pk_B_ri[j];
                ctx_CCA.b[j] -= (1 - s) * pk_B_ri[j];
        }
	}

	// Send a and b from mod q to mod p
	for(int i = 0; i < LWE_N; ++i) { ctx_CCA.a[i] = (ctx_CCA.a[i] + RD_ADD) & RD_AND; }
	for(int i = 0; i < LWE_L; ++i) { ctx_CCA.b[i] = (ctx_CCA.b[i] + RD_ADD) & RD_AND; }
	finish = clock();
	finish_cycle= rdtsc();
	elapsed1 += (finish-start);
	cycles1 += (finish_cycle - start_cycle);
	free(hash);
}
	
int Dec_CCA(){
	uint16_t decomp_delta[LWE_L];	
	int res = 0;

	start = clock();
	start_cycle = rdtsc();

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


	uint64_t delta[4] = {};
	for(int i = 0; i < 4; ++i){ for(int j = 0; j < 64; ++j){
		uint64_t a = ((uint64_t) decomp_delta[64 * i + j]) << j;
		delta[i] ^= a;
	}}
	// decomp_delta[64 * i + j] = (delta[i] >> j) & 0x01;

  	uint64_t *hash = (uint64_t *)calloc(32, sizeof(uint64_t)); //Set seed to be 5 by 5 matrix of uint64_t

	memset((unsigned char*) hash, 0, 256);
	memcpy((unsigned char*) hash, (unsigned char*) delta, 32);
	shake128((unsigned char*)hash, 256, (unsigned char*)hash, 32);

	for(int i = 0; i < 4; ++i){ ptx_CCA[i] = ctx_CCA.c1[i] ^ hash[i]; }

	if((hash[4] != ctx_CCA.c3[0]) || (hash[5] != ctx_CCA.c3[1]) || (hash[6] != ctx_CCA.c3[2])|| 
		(hash[7] != ctx_CCA.c3[3])){ res =  1; }


	uint16_t c2h_b[LWE_L];
	for(int i = 0; i < LWE_L; ++i){ c2h_b[i] = decomp_delta[i] << 15; }
	uint16_t c2h_a[LWE_N] = {};

	uint16_t r_idx[HR];
	size_t neg_start = 0;
	
	for(int i = 8; i < 30; ++i){
		int j = 6 * (i - 8);
		r_idx[j] = (uint16_t) hash[i] & 0x03ff;
		r_idx[j+1] = (uint16_t) (hash[i] >> 10) & 0x03ff;
		r_idx[j+2] = (uint16_t) (hash[i] >> 20) & 0x03ff;
		r_idx[j+3] = (uint16_t) (hash[i] >> 30) & 0x03ff;
		r_idx[j+4] = (uint16_t) (hash[i] >> 40) & 0x03ff;
		r_idx[j+5] = (uint16_t) (hash[i] >> 50) & 0x03ff;
		neg_start += (size_t) ((hash[i] >> 60) & 0x01);
		neg_start += (size_t) ((hash[i] >> 61) & 0x01);
		neg_start += (size_t) ((hash[i] >> 62) & 0x01);
		neg_start += (size_t) ((hash[i] >> 63) & 0x01);
	}
	r_idx[132] = (uint16_t) hash[30] & 0x03ff;
	r_idx[133] = (uint16_t) (hash[30] >> 10) & 0x03ff;


	for(int i = 0; i < 46; ++i){
		neg_start += (size_t) ((hash[31] >> i) & 0x01);
	}

	//cout << "Hash end" << endl;
	// Compute A^T * r and B^T * r, and then add to ct.a and ct.b, respectively.
	for (size_t i = 0; i < HR; ++i) {
        uint16_t  s = (i < neg_start)?1:0;
        uint16_t* pk_A_ri = pk_CPA.A + LWE_N * r_idx[i];
        uint16_t* pk_B_ri = pk_CPA.B + LWE_L * r_idx[i];
        for (int j = 0; j < LWE_N; ++j) { 
                c2h_a[j] += s * pk_A_ri[j];
                c2h_a[j] -= (1 - s) * pk_A_ri[j];
        }
        for (int j = 0; j < LWE_L; ++j) { 
                c2h_b[j] += s * pk_B_ri[j];
                c2h_b[j] -= (1 - s) * pk_B_ri[j];
        }
	}
	
	// Send a and b from mod q to mod p
	for(int i = 0; i < LWE_N; ++i) { c2h_a[i] = (c2h_a[i] + RD_ADD) & RD_AND; }
	for(int i = 0; i < LWE_L; ++i) { c2h_b[i] = (c2h_b[i] + RD_ADD) & RD_AND; }

	for(int i = 0; i < LWE_N; ++i){ if(c2h_a[i] != ctx_CCA.a[i]) res =  2; }
	for(int i = 0; i < LWE_L; ++i){ if(c2h_b[i] != ctx_CCA.b[i]) res =  2; }

	if(res != 0) { for(int i = 0; i < 4; ++i){ ptx_CCA[i] = 0; }}
	
	finish = clock();
	finish_cycle = rdtsc();
	elapsed2 += (finish-start);
	cycles2 += (finish_cycle - start_cycle);

	free(hash);
	return res;
}


////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////   TEST FUNCTIONS   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void EncDecTest_CPA(){
    int i;

    elapsed1 = 0;
    elapsed2 = 0;
    cycles1 = 0;
    cycles2 = 0;

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
    printf("    Enc cycles: %llu\n", cycles1/testnum/iter);
    printf("    Dec cycles: %llu\n", cycles2/testnum/iter);
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
    cycles1 = 0;
   	cycles2 = 0;

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
    printf("    Enc cycles: %llu\n", cycles1/testnum/iter);
    printf("    Dec cycles: %llu\n", cycles2/testnum/iter);
}


#ifdef PARAMS_Classical
#define KEYGEN Keygen_CPA
#define ENCDECTEST EncDecTest_CPA
#endif

#ifdef PARAMS_Recommended
#define KEYGEN Keygen_CPA
#define ENCDECTEST EncDecTest_CPA
#endif

#ifdef PARAMS_Homadd
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

////////////////////////////////////////////////////////////////////////////////////////////

int main(){
	srand((unsigned)time(NULL)+(unsigned)getpid());
	printf("\n  //////////////////////////////////////////////////////////////////\n\n");
	printf("\t\t\t"PARAMNAME" Parameter\n\n");
	printf("    LWE dimension: %d, \t\tLWR dimension: %d\n", LWE_N, LWE_M);
	printf("    Plaintext dimension: %d, \t\tPlaintext Modulus: %d bits\t\n", LWE_L, LOG_T);
	printf("    Public Key modulus: %d bits, \tCiphertext modulus: %d bits\t\n\n", LOG_Q, LOG_P);
	printf("  //////////////////////////////////////////////////////////////////\n\n");
	printf("\t\t\tPerformance Test\n\n");

	// Key Generation
	KEYGEN();

	// Encryption and Decryption Test
	ENCDECTEST();	

	printf("\n  //////////////////////////////////////////////////////////////////\n\n");

}
