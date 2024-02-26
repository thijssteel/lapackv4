#ifndef LAPACK_C_TRSV_HPP
#define LAPACK_C_TRSV_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void lapack_c_strsv(char uplo,
                    char trans,
                    char diag,
                    lapack_idx n,
                    const float* A,
                    lapack_idx lda,
                    float* x,
                    lapack_idx incx);

void lapack_c_dtrsv(char uplo,
                    char trans,
                    char diag,
                    lapack_idx n,
                    const double* A,
                    lapack_idx lda,
                    double* x,
                    lapack_idx incx);

void lapack_c_ctrsv(char uplo,
                    char trans,
                    char diag,
                    lapack_idx n,
                    const lapack_float_complex* A,
                    lapack_idx lda,
                    lapack_float_complex* x,
                    lapack_idx incx);

void lapack_c_ztrsv(char uplo,
                    char trans,
                    char diag,
                    lapack_idx n,
                    const lapack_double_complex* A,
                    lapack_idx lda,
                    lapack_double_complex* x,
                    lapack_idx incx);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif