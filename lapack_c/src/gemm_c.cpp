#include "lapack_c/gemm.h"

#ifndef USE_FORTRAN_BLAS

    #include "lapack_cpp/base.hpp"
    #include "lapack_cpp/blas/gemm.hpp"

using namespace lapack_cpp;

    #define INSTANTIATE_C_GEMM(T, name_lower, name_upper)                  \
        void lapack_c_##name_lower(                                        \
            char transa, char transb, lapack_idx m, lapack_idx n,          \
            lapack_idx k, T alpha, const T* A, lapack_idx lda, const T* B, \
            lapack_idx ldb, T beta, T* C, lapack_idx ldc)                  \
        {                                                                  \
            Op transa_ = char2op(transa);                                  \
            Op transb_ = char2op(transb);                                  \
            ConstMatrix<T, Layout::ColMajor, lapack_idx> A_(               \
                transa_ == Op::NoTrans ? m : k,                            \
                transa_ == Op::NoTrans ? k : m, A, lda);                   \
            ConstMatrix<T, Layout::ColMajor, lapack_idx> B_(               \
                transb_ == Op::NoTrans ? k : n,                            \
                transb_ == Op::NoTrans ? n : k, B, ldb);                   \
            Matrix<T, Layout::ColMajor, lapack_idx> C_(m, n, C, ldc);      \
            lapack_cpp::gemm(transa_, transb_, alpha, A_, B_, beta, C_);   \
        }

#else  // USE_FORTRAN_BLAS

    // Define interface to Fortran BLAS
    #include "lapack_c/mangling.h"

    // Macro to define both the Fortran interface and its C wrapper
    #define INSTANTIATE_C_GEMM(T, name_lower, name_upper)                      \
        void FORTRAN_TO_C_NAME(name_lower, name_upper)(                        \
            const char* transa, const char* transb, const lapack_idx* m,       \
            const lapack_idx* n, const lapack_idx* k, const T* alpha,          \
            const T* A, const lapack_idx* lda, const T* B,                     \
            const lapack_idx* ldb, const T* beta, T* C,                        \
            const lapack_idx* ldc);                                            \
        void lapack_c_##name_lower(                                            \
            char transa, char transb, lapack_idx m, lapack_idx n,              \
            lapack_idx k, T alpha, const T* A, lapack_idx lda, const T* B,     \
            lapack_idx ldb, T beta, T* C, lapack_idx ldc)                      \
        {                                                                      \
            FORTRAN_TO_C_NAME(name_lower, name_upper)                          \
            (&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, \
             &ldc);                                                            \
        }

#endif  // USE_FORTRAN_BLAS

// Actual definitions
extern "C" {

INSTANTIATE_C_GEMM(float, sgemm, SGEMM)
INSTANTIATE_C_GEMM(double, dgemm, DGEMM)
INSTANTIATE_C_GEMM(lapack_float_complex, cgemm, CGEMM)
INSTANTIATE_C_GEMM(lapack_double_complex, zgemm, ZGEMM)

}  // extern "C"

#undef INSTANTIATE_C_GEMM