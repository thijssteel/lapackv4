#include "lapack_c/gemm.h"

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/blas/gemm.hpp"

using namespace lapack_cpp;

namespace lapack_c {

/**
 * Wrapper for the C++ gemm function
 */
template <typename T, typename idx_t>
void gemm_cpp_wrapper(char transa,
                      char transb,
                      idx_t m,
                      idx_t n,
                      idx_t k,
                      T alpha,
                      const T* A,
                      idx_t lda,
                      const T* B,
                      idx_t ldb,
                      T beta,
                      T* C,
                      idx_t ldc)
{
    // Construct wrapper objects
    Op transa_ = char2op(transa);
    Op transb_ = char2op(transb);
    ConstMatrix<T, idx_t> A_(transa_==Op::NoTrans ? m : k, transa_==Op::NoTrans ? k : m, A, lda);
    ConstMatrix<T, idx_t> B_(transb_==Op::NoTrans ? k : n, transb_==Op::NoTrans ? n : k, B, ldb);
    Matrix<T, idx_t> C_(m, n, C, ldc);

    // Call the C++ function
    lapack_cpp::gemm( transa_, transb_, alpha, A_, B_, beta, C_);
}

extern "C" {
//
// Define actual wrappers
//
void lapack_c_dgemm(char transa,
                    char transb,
                    int m,
                    int n,
                    int k,
                    double alpha,
                    const double* A,
                    int lda,
                    const double* B,
                    int ldb,
                    double beta,
                    double* C,
                    int ldc)
{
    gemm_cpp_wrapper(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

void lapack_c_sgemm(char transa,
                    char transb,
                    int m,
                    int n,
                    int k,
                    float alpha,
                    const float* A,
                    int lda,
                    const float* B,
                    int ldb,
                    float beta,
                    float* C,
                    int ldc)
{
    gemm_cpp_wrapper(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

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
                    lapack_idx ldc)
{
    gemm_cpp_wrapper(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

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
                    lapack_idx ldc)
{
    gemm_cpp_wrapper(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

}  // extern "C"

}  // namespace lapack_c