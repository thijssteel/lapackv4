#ifndef LAPACK_C_GEMM_HPP
#define LAPACK_C_GEMM_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void lapack_c_dgemm(char transa,
                    char transb,
                    lapack_idx m,
                    lapack_idx n,
                    lapack_idx k,
                    double alpha,
                    const double* A,
                    lapack_idx lda,
                    const double* B,
                    lapack_idx ldb,
                    double beta,
                    double* C,
                    lapack_idx ldc);

void lapack_c_sgemm(char transa,
                    char transb,
                    lapack_idx m,
                    lapack_idx n,
                    lapack_idx k,
                    float alpha,
                    const float* A,
                    lapack_idx lda,
                    const float* B,
                    lapack_idx ldb,
                    float beta,
                    float* C,
                    lapack_idx ldc);

void lapack_c_cgemm(char transa,
                    char transb,
                    lapack_idx m,
                    lapack_idx n,
                    lapack_idx k,
                    lapack_float_complex alpha,
                    const lapack_float_complex* A,
                    lapack_idx lda,
                    const lapack_float_complex* B,
                    lapack_idx ldb,
                    lapack_float_complex beta,
                    lapack_float_complex* C,
                    lapack_idx ldc);

void lapack_c_zgemm(char transa,
                    char transb,
                    lapack_idx m,
                    lapack_idx n,
                    lapack_idx k,
                    lapack_double_complex alpha,
                    const lapack_double_complex* A,
                    lapack_idx lda,
                    const lapack_double_complex* B,
                    lapack_idx ldb,
                    lapack_double_complex beta,
                    lapack_double_complex* C,
                    lapack_idx ldc);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif