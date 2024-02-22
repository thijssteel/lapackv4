#ifndef LAPACK_CPP_GEMM_HPP
#define LAPACK_CPP_GEMM_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"
namespace lapack_cpp {

/**
 * Perform a matrix-matrix multiplication.
 *
 * @tparam T
 * @param A
 * @param B
 * @param C
 */
template <typename T, typename idx_t, Layout layout>
void gemm(Op transA,
          Op transB,
          T alpha,
          const ConstMatrix<T, idx_t, layout>& A,
          const ConstMatrix<T, idx_t, layout>& B,
          T beta,
          const Matrix<T, idx_t, layout>& C)
{
    // TODO: Add checks for the dimensions of A, B, and C
    T zero = 0.;
    T one = 1.;

    const idx_t m = (transA == Op::NoTrans) ? A.num_rows() : A.num_columns();
    const idx_t n = (transB == Op::NoTrans) ? B.num_columns() : B.num_rows();
    const idx_t k = (transA == Op::NoTrans) ? A.num_columns() : A.num_rows();

    if (m == 0 || n == 0 || k == 0) return;

    if (alpha == zero) {
        if (beta == zero) {
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i)
                    C(i, j) = zero;
            }
        }
        else if (beta != one) {
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i)
                    C(i, j) *= beta;
            }
        }
        return;
    }

    // alpha != zero
    if (transA == Op::NoTrans) {
        if (transB == Op::NoTrans) {
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i)
                    C(i, j) *= beta;
                for (idx_t l = 0; l < k; ++l) {
                    T alpha_Blj = alpha * B(l, j);
                    for (idx_t i = 0; i < m; ++i)
                        C(i, j) += A(i, l) * alpha_Blj;
                }
            }
        }
        else if (transB == Op::Trans) {
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i)
                    C(i, j) *= beta;
                for (idx_t l = 0; l < k; ++l) {
                    T alpha_Bjl = alpha * B(j, l);
                    for (idx_t i = 0; i < m; ++i)
                        C(i, j) += A(i, l) * alpha_Bjl;
                }
            }
        }
        else {  // transB == Op::ConjTrans
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i)
                    C(i, j) *= beta;
                for (idx_t l = 0; l < k; ++l) {
                    T alpha_Bjl = alpha * conj(B(j, l));
                    for (idx_t i = 0; i < m; ++i)
                        C(i, j) += A(i, l) * alpha_Bjl;
                }
            }
        }
    }
    else if (transA == Op::Trans) {
        if (transB == Op::NoTrans) {
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i) {
                    T sum = zero;
                    for (idx_t l = 0; l < k; ++l)
                        sum += A(l, i) * B(l, j);
                    C(i, j) = alpha * sum + beta * C(i, j);
                }
            }
        }
        else if (transB == Op::Trans) {
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i) {
                    T sum = zero;
                    for (idx_t l = 0; l < k; ++l)
                        sum += A(l, i) * B(j, l);
                    C(i, j) = alpha * sum + beta * C(i, j);
                }
            }
        }
        else {  // transB == Op::ConjTrans
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i) {
                    T sum = zero;
                    for (idx_t l = 0; l < k; ++l)
                        sum += A(l, i) * conj(B(j, l));
                    C(i, j) = alpha * sum + beta * C(i, j);
                }
            }
        }
    }
    else {  // transA == Op::ConjTrans
        if (transB == Op::NoTrans) {
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i) {
                    T sum = zero;
                    for (idx_t l = 0; l < k; ++l)
                        sum += conj(A(l, i)) * B(l, j);
                    C(i, j) = alpha * sum + beta * C(i, j);
                }
            }
        }
        else if (transB == Op::Trans) {
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i) {
                    T sum = zero;
                    for (idx_t l = 0; l < k; ++l)
                        sum += conj(A(l, i)) * B(j, l);
                    C(i, j) = alpha * sum + beta * C(i, j);
                }
            }
        }
        else {  // transB == Op::ConjTrans
            for (idx_t j = 0; j < n; ++j) {
                for (idx_t i = 0; i < m; ++i) {
                    T sum = zero;
                    for (idx_t l = 0; l < k; ++l)
                        sum += A(l, i) * B(j, l);
                    C(i, j) = alpha * conj(sum) + beta * C(i, j);
                }
            }
        }
    }
}

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_GEMM_HPP