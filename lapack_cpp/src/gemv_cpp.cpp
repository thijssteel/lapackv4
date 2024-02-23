#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/blas/gemv.hpp"

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
inline void gemv_c_wrapper(Op trans,
                           T alpha,
                           const ConstMatrix<T, idx_t, layout>& A,
                           const ConstVector<T, idx_t>& x,
                           T beta,
                           const Vector<T, idx_t>& y)
{
    static_assert(layout == Layout::ColMajor);

    assert(((trans == Op::NoTrans) ? A.num_columns() : A.num_rows()) ==
           x.size());
    assert(((trans == Op::NoTrans) ? A.num_rows() : A.num_columns()) ==
           y.size());

    if constexpr (std::is_same<T, double>::value) {
        lapack_c_dgemv((char)trans, A.num_rows(), A.num_columns(), alpha,
                       A.ptr(), A.ldim(), x.ptr(), x.stride(), beta, y.ptr(),
                       y.stride());
    }
    else if constexpr (std::is_same<T, float>::value) {
        lapack_c_sgemv((char)trans, A.num_rows(), A.num_columns(), alpha,
                       A.ptr(), A.ldim(), x.ptr(), x.stride(), beta, y.ptr(),
                       y.stride());
    }
    else if constexpr (std::is_same<T, std::complex<float>>::value) {
        lapack_c_cgemv((char)trans, A.num_rows(), A.num_columns(), alpha,
                       A.ptr(), A.ldim(), x.ptr(), x.stride(), beta, y.ptr(),
                       y.stride());
    }
    else if constexpr (std::is_same<T, std::complex<double>>::value) {
        lapack_c_zgemv((char)trans, A.num_rows(), A.num_columns(), alpha,
                       A.ptr(), A.ldim(), x.ptr(), x.stride(), beta, y.ptr(),
                       y.stride());
    }
    else {
        assert(false);
    }
}

// We have a lot of types to instantiate for, so we use a macro to avoid
// repetition.
#define INSTANTIATE_GEMV(T, idx_t, layout)                               \
    template <>                                                          \
    void gemv(Op trans, T alpha, const ConstMatrix<T, idx_t, layout>& A, \
              const ConstVector<T, idx_t>& x, T beta,                    \
              const Vector<T, idx_t>& y)                                 \
    {                                                                    \
        gemv_c_wrapper(trans, alpha, A, x, beta, y);                     \
    }

INSTANTIATE_GEMV(float, size_t, Layout::ColMajor)
INSTANTIATE_GEMV(double, size_t, Layout::ColMajor)
INSTANTIATE_GEMV(std::complex<float>, size_t, Layout::ColMajor)
INSTANTIATE_GEMV(std::complex<double>, size_t, Layout::ColMajor)
INSTANTIATE_GEMV(float, int, Layout::ColMajor)
INSTANTIATE_GEMV(double, int, Layout::ColMajor)
INSTANTIATE_GEMV(std::complex<float>, int, Layout::ColMajor)
INSTANTIATE_GEMV(std::complex<double>, int, Layout::ColMajor)

// Row major is not yet supported so don't instantiate it

#undef INSTANTIATE_GEMV

}  // namespace lapack_cpp