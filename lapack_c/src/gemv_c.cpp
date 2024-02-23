#include "lapack_c/gemv.h"

// Note, these is no C++ implementation for GEMV yet, so always use the Fortran
// interface

#include <iostream>

#include "lapack_c/mangling.h"

// Macro to define both the Fortran interface and its C wrapper
#define INSTANTIATE_C_GEMV(T, name_lower, name_upper)                         \
    void FORTRAN_TO_C_NAME(name_lower, name_upper)(                           \
        const char* trans, const lapack_idx* m, const lapack_idx* n,          \
        const T* alpha, const T* A, const lapack_idx* lda, const T* x,        \
        const lapack_idx* incx, const T* beta, T* y, const lapack_idx* incy); \
    void lapack_c_##name_lower(char trans, lapack_idx m, lapack_idx n,        \
                               T alpha, const T* A, lapack_idx lda,           \
                               const T* x, lapack_idx incx, T beta, T* y,     \
                               lapack_idx incy)                               \
    {                                                                         \
        FORTRAN_TO_C_NAME(name_lower, name_upper)                             \
        (&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);         \
    }

// Actual definitions
extern "C" {

INSTANTIATE_C_GEMV(float, sgemv, SGEMV)
INSTANTIATE_C_GEMV(double, dgemv, DGEMV)
INSTANTIATE_C_GEMV(lapack_float_complex, cgemv, CGEMV)
INSTANTIATE_C_GEMV(lapack_double_complex, zgemv, ZGEMV)

}  // extern "C"

#undef INSTANTIATE_C_GEMV