#include <string.h>
#include <assert.h>

#include "../CTIDH/ctidh.h"
#include "primes.h"
#include "mont.h"
#include "elligator.h"
// #include "fp2.h"

#ifdef ENABLE_CT_TESTING
#include <valgrind/memcheck.h>
#endif

// static void fp2_print(fp2 a){
//     printf("0x%016lX,", A[i]);
    
//     fp_print(a.re);
//     fp_print(a.im);
// }


// For computing [(p + 1) / l_i]P, i:=0, ..., (N - 1)
void cofactor_multiples(proj P[], proj const A, size_t lower, size_t upper)
{
    assert(lower < upper);
    if (upper - lower == 1)
        return;

    int i;
    size_t mid = lower + (upper - lower + 1) / 2;

    // proj_copy(P[mid], (const fp*)P[lower]);
    fp_copy(P[mid].x, P[lower].x);
    fp_copy(P[mid].z, P[lower].z);

    for (i = lower; i < (int)mid; i++)
        xMUL_dac(&P[mid], &A, 1, &P[mid], primes_dac[i], primes_daclen[i], primes_daclen[i]);
    // xmul(P[mid], i, (const fp*)P[mid], A);

    for (i = (int)mid; i < (int)upper; i++)
        xMUL_dac(&P[lower], &A, 1, &P[lower], primes_dac[i], primes_daclen[i], primes_daclen[i]);
    // xmul(P[lower], i, (const fp*)P[lower], A);

    cofactor_multiples(P, A, lower, mid);
    cofactor_multiples(P, A, mid, upper);
}


static void clearpublicprimes_vali(proj *P, const proj *A24)
{
  // clear powers of 2
  for (int64_t i = 0; i < two_cofactor; i++)
  {
    xDBL(P, P, A24, 0);
  }
}

// wombat validate by checking full order point
bool validate(public_key const *in){
    proj A, A24;
    fp_copy(A.x, in->A);
    fp_copy(A.z, fp_1);
    xA24(&A24, &A);

    proj Tp, Tm, Pp[primes_num], Pm[primes_num];

    uint8_t boolp = 0, boolm = 0;

#ifdef ENABLE_CT_TESTING
    VALGRIND_MAKE_MEM_DEFINED(&A, sizeof(proj));
    VALGRIND_MAKE_MEM_DEFINED(&boolp, sizeof(uint8_t));
    VALGRIND_MAKE_MEM_DEFINED(&boolm, sizeof(uint8_t));
    VALGRIND_MAKE_MEM_DEFINED(&Tp, sizeof(proj));
    VALGRIND_MAKE_MEM_DEFINED(Pp, sizeof(proj) * primes_num);
    VALGRIND_MAKE_MEM_DEFINED(Pm, sizeof(proj) * primes_num);
#endif
    fp seed;
    fp_set(seed, in->seed);
    fp_enc(seed, seed);
    elligator_seeded(&Tp, &Tm, &A, (const fp *)seed);

    clearpublicprimes_vali(&Tp, &A);
    clearpublicprimes_vali(&Tm, &A);

    int64_t primes_needed = batch_stop[primes_batches - 1] + 1;

    int64_t j;
    // clear ells not used in the keyspace
    for (j = primes_needed; j < primes_num; j++)
    {
        xMUL_dac(&Tp, &A, 0, &Tp, primes_dac[j], primes_daclen[j], primes_daclen[j]);
        xMUL_dac(&Tm, &A, 0, &Tm, primes_dac[j], primes_daclen[j], primes_daclen[j]);
    }

    // Checking if Tp is an order (p+1)/(2^e)
    proj_copy(&Pp[0], &Tp);
    cofactor_multiples(Pp, A, 0, primes_needed);
    boolp = 1;
    for (j = batch_start[0]; j < primes_needed; j++)
        boolp &= (1 - fp_iszero(Pp[j].z));


    proj_copy(&Pm[0], &Tm);
    cofactor_multiples(Pm, A, 0, primes_needed);

    boolm = 1;
    for (j = batch_start[0]; j < primes_needed; j++)
        boolm &= (1 - fp_iszero(Pm[j].z));
    
    return boolp & boolm;
}
