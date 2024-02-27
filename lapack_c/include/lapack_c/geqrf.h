#ifndef LAPACK_C_GEQRF_HPP
#define LAPACK_C_GEQRF_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void lapack_c_sgeqrf(lapack_idx m,
                     lapack_idx n,
                     float* A,
                     lapack_idx lda,
                     float* tau,
                     float* work,
                     lapack_idx lwork,
                     lapack_idx* info);

void lapack_c_dgeqrf(lapack_idx m,
                     lapack_idx n,
                     double* A,
                     lapack_idx lda,
                     double* tau,
                     double* work,
                     lapack_idx lwork,
                     lapack_idx* info);

void lapack_c_cgeqrf(lapack_idx m,
                     lapack_idx n,
                     lapack_float_complex* A,
                     lapack_idx lda,
                     lapack_float_complex* tau,
                     lapack_float_complex* work,
                     lapack_idx lwork,
                     lapack_idx* info);

void lapack_c_zgeqrf(lapack_idx m,
                     lapack_idx n,
                     lapack_double_complex* A,
                     lapack_idx lda,
                     lapack_double_complex* tau,
                     lapack_double_complex* work,
                     lapack_idx lwork,
                     lapack_idx* info);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif