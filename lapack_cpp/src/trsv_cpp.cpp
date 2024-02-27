#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/blas/trsv.hpp"

namespace lapack_cpp {

/**
 * Templated wrapper around c functions, will be instantiated for each type.
 *
 * @tparam T
 * @param A
 * @param B
 * @param C
 */
template <typename T, Layout layout, typename idx_t>
inline void trsv_c_wrapper(Uplo uplo,
                           Op trans,
                           Diag diag,
                           const ConstMatrix<T, layout, idx_t>& A,
                           const Vector<T, idx_t>& x)
{
    static_assert(layout == Layout::ColMajor);

    assert(A.num_rows() == A.num_columns());
    assert(A.num_rows() == x.size());

    if constexpr (std::is_same<T, double>::value) {
        lapack_c_dtrsv((char)uplo, (char)trans, (char)diag, A.num_columns(), A.ptr(),
                       A.ldim(), x.ptr(), x.stride());
    }
    else if constexpr (std::is_same<T, float>::value) {
        lapack_c_strsv((char)uplo, (char)trans, (char)diag, A.num_columns(), A.ptr(),
                       A.ldim(), x.ptr(), x.stride());
    }
    else if constexpr (std::is_same<T, std::complex<float>>::value) {
        lapack_c_ctrsv((char)uplo, (char)trans, (char)diag, A.num_columns(), A.ptr(),
                       A.ldim(), x.ptr(), x.stride());
    }
    else if constexpr (std::is_same<T, std::complex<double>>::value) {
        lapack_c_ztrsv((char)uplo, (char)trans, (char)diag, A.num_columns(), A.ptr(),
                       A.ldim(), x.ptr(), x.stride());
    }
    else {
        assert(false);
    }
}

// We have a lot of types to instantiate for, so we use a macro to avoid
// repetition.
#define INSTANTIATE_TRSV(T, idx_t, layout)            \
    template <>                                       \
    void trsv(Uplo uplo, Op trans, Diag diag,         \
              const ConstMatrix<T, layout, idx_t>& A, \
              const Vector<T, idx_t>& x)              \
    {                                                 \
        trsv_c_wrapper(uplo, trans, diag, A, x);  \
    }

INSTANTIATE_TRSV(float, lapack_idx_t, Layout::ColMajor)
INSTANTIATE_TRSV(double, lapack_idx_t, Layout::ColMajor)
INSTANTIATE_TRSV(std::complex<float>, lapack_idx_t, Layout::ColMajor)
INSTANTIATE_TRSV(std::complex<double>, lapack_idx_t, Layout::ColMajor)
INSTANTIATE_TRSV(float, int, Layout::ColMajor)
INSTANTIATE_TRSV(double, int, Layout::ColMajor)
INSTANTIATE_TRSV(std::complex<float>, int, Layout::ColMajor)
INSTANTIATE_TRSV(std::complex<double>, int, Layout::ColMajor)

// Row major is not yet supported so don't instantiate it

#undef INSTANTIATE_TRSV

}  // namespace lapack_cpp