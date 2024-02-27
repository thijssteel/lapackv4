#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/lapack/geqrf.hpp"

namespace lapack_cpp {

/**
 * Templated wrapper around c functions, will be instantiated for each type.
 *
 * @tparam T
 * @param A
 * @param B
 * @param C
 */
template <typename T, Layout layout, typename idx_t, bool aligned>
inline void geqrf_c_wrapper(const Matrix<T, layout, idx_t>& A,
                            const Vector<T, idx_t>& tau,
                            const MemoryBlock<T, idx_t, aligned>& work)
{
    static_assert(layout == Layout::ColMajor);

    int info;

    if constexpr (std::is_same<T, double>::value) {
        lapack_c_dgeqrf(A.num_rows(), A.num_columns(), A.ptr(), A.ldim(),
                        tau.ptr(), work.ptr(), work.size(), &info);
    }
    else if constexpr (std::is_same<T, float>::value) {
        lapack_c_sgeqrf(A.num_rows(), A.num_columns(), A.ptr(), A.ldim(),
                        tau.ptr(), work.ptr(), work.size(), &info);
    }
    else if constexpr (std::is_same<T, std::complex<float>>::value) {
        lapack_c_cgeqrf(A.num_rows(), A.num_columns(), A.ptr(), A.ldim(),
                        tau.ptr(), work.ptr(), work.size(), &info);
    }
    else if constexpr (std::is_same<T, std::complex<double>>::value) {
        lapack_c_zgeqrf(A.num_rows(), A.num_columns(), A.ptr(), A.ldim(),
                        tau.ptr(), work.ptr(), work.size(), &info);
    }
    else {
        assert(false);
    }

    // TODO: check for info and throw exception
}

/**
 * Templated wrapper around c functions, will be instantiated for each type.
 *
 * @tparam T
 * @param A
 * @param B
 * @param C
 */
template <typename T, Layout layout, typename idx_t>
inline idx_t geqrf_work_c_wrapper(const Matrix<T, layout, idx_t>& A,
                                  const Vector<T, idx_t>& tau)
{
    static_assert(layout == Layout::ColMajor);

    int info;
    T dummywork;

    if constexpr (std::is_same<T, double>::value) {
        lapack_c_dgeqrf(A.num_rows(), A.num_columns(), A.ptr(), A.ldim(),
                        tau.ptr(), &dummywork, -1, &info);
    }
    else if constexpr (std::is_same<T, float>::value) {
        lapack_c_sgeqrf(A.num_rows(), A.num_columns(), A.ptr(), A.ldim(),
                        tau.ptr(), &dummywork, -1, &info);
    }
    else if constexpr (std::is_same<T, std::complex<float>>::value) {
        lapack_c_cgeqrf(A.num_rows(), A.num_columns(), A.ptr(), A.ldim(),
                        tau.ptr(), &dummywork, -1, &info);
    }
    else if constexpr (std::is_same<T, std::complex<double>>::value) {
        lapack_c_zgeqrf(A.num_rows(), A.num_columns(), A.ptr(), A.ldim(),
                        tau.ptr(), &dummywork, -1, &info);
    }
    else {
        assert(false);
    }

    // TODO: check for info and throw exception

    return static_cast<idx_t>(real(dummywork));
}

// We have a lot of types to instantiate for, so we use a macro to avoid
// repetition.
#define INSTANTIATE_GEQRF(T, idx_t, layout, aligned)                           \
    template <>                                                                \
    void geqrf(const Matrix<T, layout, idx_t>& A, const Vector<T, idx_t>& tau, \
               const MemoryBlock<T, idx_t, aligned>& work)                     \
    {                                                                          \
        geqrf_c_wrapper(A, tau, work);                                         \
    }

INSTANTIATE_GEQRF(float, lapack_idx_t, Layout::ColMajor, true)
INSTANTIATE_GEQRF(double, lapack_idx_t, Layout::ColMajor, true)
INSTANTIATE_GEQRF(std::complex<float>, lapack_idx_t, Layout::ColMajor, true)
INSTANTIATE_GEQRF(std::complex<double>, lapack_idx_t, Layout::ColMajor, true)
INSTANTIATE_GEQRF(float, int, Layout::ColMajor, true)
INSTANTIATE_GEQRF(double, int, Layout::ColMajor, true)
INSTANTIATE_GEQRF(std::complex<float>, int, Layout::ColMajor, true)
INSTANTIATE_GEQRF(std::complex<double>, int, Layout::ColMajor, true)

INSTANTIATE_GEQRF(float, lapack_idx_t, Layout::ColMajor, false)
INSTANTIATE_GEQRF(double, lapack_idx_t, Layout::ColMajor, false)
INSTANTIATE_GEQRF(std::complex<float>, lapack_idx_t, Layout::ColMajor, false)
INSTANTIATE_GEQRF(std::complex<double>, lapack_idx_t, Layout::ColMajor, false)
INSTANTIATE_GEQRF(float, int, Layout::ColMajor, false)
INSTANTIATE_GEQRF(double, int, Layout::ColMajor, false)
INSTANTIATE_GEQRF(std::complex<float>, int, Layout::ColMajor, false)
INSTANTIATE_GEQRF(std::complex<double>, int, Layout::ColMajor, false)

// Row major is not yet supported so don't instantiate it

#undef INSTANTIATE_GEQRF

#define INSTANTIATE_GEQRF_WORKQUERY(T, idx_t, layout)        \
    template <>                                              \
    idx_t geqrf_workquery(const Matrix<T, layout, idx_t>& A, \
                          const Vector<T, idx_t>& tau)       \
    {                                                        \
        return geqrf_work_c_wrapper(A, tau);                 \
    }

INSTANTIATE_GEQRF_WORKQUERY(float, lapack_idx_t, Layout::ColMajor)
INSTANTIATE_GEQRF_WORKQUERY(double, lapack_idx_t, Layout::ColMajor)
INSTANTIATE_GEQRF_WORKQUERY(std::complex<float>, lapack_idx_t, Layout::ColMajor)
INSTANTIATE_GEQRF_WORKQUERY(std::complex<double>,
                            lapack_idx_t,
                            Layout::ColMajor)
INSTANTIATE_GEQRF_WORKQUERY(float, int, Layout::ColMajor)
INSTANTIATE_GEQRF_WORKQUERY(double, int, Layout::ColMajor)
INSTANTIATE_GEQRF_WORKQUERY(std::complex<float>, int, Layout::ColMajor)
INSTANTIATE_GEQRF_WORKQUERY(std::complex<double>, int, Layout::ColMajor)

#undef INSTANTIATE_GEQRF_WORKQUERY

}  // namespace lapack_cpp