#ifdef USE_FORTRAN_BLAS

#include "lapack_cpp/blas/gemm.hpp"

#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"

namespace lapack_cpp {

/**
 * Templated wrapper around c functions, will be instantiated for each type.
 *
 * @tparam T
 * @param A
 * @param B
 * @param C
 */
template <typename T, typename idx_t, Layout layout>
inline void gemm_c_wrapper(Op transA,
                           Op transB,
                           T alpha,
                           const ConstMatrix<T, idx_t, layout>& A,
                           const ConstMatrix<T, idx_t, layout>& B,
                           T beta,
                           const Matrix<T, idx_t, layout>& C)
{
    if (layout == Layout::RowMajor) {
        // Row-major -> Col-major = Transpose
        // C = alpha * A * B + beta * C -> C^T = alpha * B^T * A^T + beta * C^T
        // Reinterpret A and B as column-major and call gemm with switched
        // arguments
        ConstMatrix<T, idx_t, Layout::ColMajor> A_t(
            A.num_columns(), A.num_rows(), A.ptr(), A.ldim());
        ConstMatrix<T, idx_t, Layout::ColMajor> B_t(
            B.num_columns(), B.num_rows(), B.ptr(), B.ldim());
        Matrix<T, idx_t, Layout::ColMajor> C_t(C.num_columns(), C.num_rows(),
                                               C.ptr(), C.ldim());
        gemm_c_wrapper(transB, transA, alpha, B_t, A_t, beta, C_t);
        return;
    }

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
        lapack_c_dgemm((char)transA, (char)transB, m, n, k, alpha, A.ptr(),
                       A.ldim(), B.ptr(), B.ldim(), beta, C.ptr(), C.ldim());
    }
    else if constexpr (std::is_same<T, float>::value) {
        lapack_c_sgemm((char)transA, (char)transB, m, n, k, alpha, A.ptr(),
                       A.ldim(), B.ptr(), B.ldim(), beta, C.ptr(), C.ldim());
    }
    else if constexpr (std::is_same<T, std::complex<float>>::value) {
        lapack_c_cgemm((char)transA, (char)transB, m, n, k, alpha, A.ptr(),
                       A.ldim(), B.ptr(), B.ldim(), beta, C.ptr(), C.ldim());
    }
    else if constexpr (std::is_same<T, std::complex<double>>::value) {
        lapack_c_zgemm((char)transA, (char)transB, m, n, k, alpha, A.ptr(),
                       A.ldim(), B.ptr(), B.ldim(), beta, C.ptr(), C.ldim());
    }
    else {
        assert(false);
    }
}


// We have a lot of types to instantiate for, so we use a macro to avoid
// repetition.
#define INSTANTIATE_GEMM(T, idx_t, layout)                    \
    template <>                                               \
    void gemm(Op transA, Op transB, T alpha,                  \
              const ConstMatrix<T, idx_t, layout>& A,         \
              const ConstMatrix<T, idx_t, layout>& B, T beta, \
              const Matrix<T, idx_t, layout>& C)              \
    {                                                         \
        gemm_c_wrapper(transA, transB, alpha, A, B, beta, C); \
    }

INSTANTIATE_GEMM(float, size_t, Layout::ColMajor)
INSTANTIATE_GEMM(double, size_t, Layout::ColMajor)
INSTANTIATE_GEMM(std::complex<float>, size_t, Layout::ColMajor)
INSTANTIATE_GEMM(std::complex<double>, size_t, Layout::ColMajor)
INSTANTIATE_GEMM(float, int, Layout::ColMajor)
INSTANTIATE_GEMM(double, int, Layout::ColMajor)
INSTANTIATE_GEMM(std::complex<float>, int, Layout::ColMajor)
INSTANTIATE_GEMM(std::complex<double>, int, Layout::ColMajor)

INSTANTIATE_GEMM(float, size_t, Layout::RowMajor)
INSTANTIATE_GEMM(double, size_t, Layout::RowMajor)
INSTANTIATE_GEMM(std::complex<float>, size_t, Layout::RowMajor)
INSTANTIATE_GEMM(std::complex<double>, size_t, Layout::RowMajor)
INSTANTIATE_GEMM(float, int, Layout::RowMajor)
INSTANTIATE_GEMM(double, int, Layout::RowMajor)
INSTANTIATE_GEMM(std::complex<float>, int, Layout::RowMajor)
INSTANTIATE_GEMM(std::complex<double>, int, Layout::RowMajor)

#undef INSTANTIATE_GEMM

}  // namespace lapack_cpp

#endif  // USE_FORTRAN_BLAS