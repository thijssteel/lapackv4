#ifndef LAPACK_C_GEMV_HPP
#define LAPACK_C_GEMV_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void lapack_c_dgemv(char trans,
                    lapack_idx m,
                    lapack_idx n,
                    double alpha,
                    const double* A,
                    lapack_idx lda,
                    const double* x,
                    lapack_idx incx,
                    double beta,
                    double* y,
                    lapack_idx incy);

void lapack_c_sgemv(char trans,
                    lapack_idx m,
                    lapack_idx n,
                    float alpha,
                    const float* A,
                    lapack_idx lda,
                    const float* x,
                    lapack_idx incx,
                    float beta,
                    float* y,
                    lapack_idx incy);

void lapack_c_cgemv(char trans,
                    lapack_idx m,
                    lapack_idx n,
                    lapack_float_complex alpha,
                    const lapack_float_complex* A,
                    lapack_idx lda,
                    const lapack_float_complex* x,
                    lapack_idx incx,
                    lapack_float_complex beta,
                    lapack_float_complex* y,
                    lapack_idx incy);

void lapack_c_zgemv(char trans,
                    lapack_idx m,
                    lapack_idx n,
                    lapack_double_complex alpha,
                    const lapack_double_complex* A,
                    lapack_idx lda,
                    const lapack_double_complex* x,
                    lapack_idx incx,
                    lapack_double_complex beta,
                    lapack_double_complex* y,
                    lapack_idx incy);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif