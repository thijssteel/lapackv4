#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/blas/gemm.hpp"

namespace lapack_cpp {

/**
 * Templated wrapper around c functions, will be instantiated for each type.
 *
 * @tparam T
 * @param A
 * @param B
 * @param C
 */
template <typename T, typename idx_t>
inline void gemm_c_wrapper(Op transA,
                           Op transB,
                           T alpha,
                           const ConstMatrix<T, idx_t>& A,
                           const ConstMatrix<T, idx_t>& B,
                           T beta,
                           const Matrix<T, idx_t>& C)
{
    assert(((transA == Op::NoTrans) ? A.num_columns() : A.num_rows()) ==
           ((transB == Op::NoTrans) ? B.num_rows() : B.num_columns()));
    assert(((transA == Op::NoTrans) ? A.num_rows() : A.num_columns()) ==
           C.num_rows());
    assert(((transB == Op::NoTrans) ? B.num_columns() : B.num_rows()) ==
           C.num_columns());

    const lapack_idx m =
        (transA == Op::NoTrans) ? A.num_rows() : A.num_columns();
    const lapack_idx n =
        (transB == Op::NoTrans) ? B.num_columns() : B.num_rows();
    const lapack_idx k =
        (transA == Op::NoTrans) ? A.num_columns() : A.num_rows();

    if constexpr (std::is_same<T, double>::value) {
        lapack_c_dgemm((char)transA, (char)transB, m, n, k, alpha, A.ptr(), A.ldim(),
                       B.ptr(), B.ldim(), beta, C.ptr(), C.ldim());
    }
    else if constexpr (std::is_same<T, float>::value) {
        lapack_c_sgemm((char)transA, (char)transB, m, n, k, alpha, A.ptr(), A.ldim(),
                       B.ptr(), B.ldim(), beta, C.ptr(), C.ldim());
    }
    else if constexpr (std::is_same<T, std::complex<float>>::value) {
        lapack_c_cgemm((char)transA, (char)transB, m, n, k, alpha, A.ptr(), A.ldim(),
                       B.ptr(), B.ldim(), beta, C.ptr(), C.ldim());
    }
    else if constexpr (std::is_same<T, std::complex<double>>::value) {
        lapack_c_zgemm((char)transA, (char)transB, m, n, k, alpha, A.ptr(), A.ldim(),
                       B.ptr(), B.ldim(), beta, C.ptr(), C.ldim());
    }
    else {
        assert(false);
    }
}

// Specialization for double, size_t
template <>
void gemm(Op transA,
          Op transB,
          double alpha,
          const ConstMatrix<double, size_t>& A,
          const ConstMatrix<double, size_t>& B,
          double beta,
          const Matrix<double, size_t>& C)
{
    gemm_c_wrapper(transA, transB, alpha, A, B, beta, C);
}

// Specialization for float, size_t
template <>
void gemm(Op transA,
          Op transB,
          float alpha,
          const ConstMatrix<float, size_t>& A,
          const ConstMatrix<float, size_t>& B,
          float beta,
          const Matrix<float, size_t>& C)
{
    gemm_c_wrapper(transA, transB, alpha, A, B, beta, C);
}

// Specialization for std::complex<float>, size_t
template <>
void gemm(Op transA,
          Op transB,
          std::complex<float> alpha,
          const ConstMatrix<std::complex<float>, size_t>& A,
          const ConstMatrix<std::complex<float>, size_t>& B,
          std::complex<float> beta,
          const Matrix<std::complex<float>, size_t>& C)
{
    gemm_c_wrapper(transA, transB, alpha, A, B, beta, C);
}

// Specialization for std::complex<double>, size_t
template <>
void gemm(Op transA,
          Op transB,
          std::complex<double> alpha,
          const ConstMatrix<std::complex<double>, size_t>& A,
          const ConstMatrix<std::complex<double>, size_t>& B,
          std::complex<double> beta,
          const Matrix<std::complex<double>, size_t>& C)
{
    gemm_c_wrapper(transA, transB, alpha, A, B, beta, C);
}

}  // namespace lapack_cpp
