#ifndef LAPACK_CPP_GEQRF_HPP
#define LAPACK_CPP_GEQRF_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"
namespace lapack_cpp {

/**
 * Compute the QR factorization of a general m-by-n matrix A.
 */
template <typename T, Layout layout, typename idx_t, bool aligned>
void geqrf(const Matrix<T, layout, idx_t>& A,
           const Vector<T, idx_t>& tau,
           const MemoryBlock<T, idx_t, aligned>& work);

/**
 * Return the optimal work array size for geqrf.
 */
template <typename T, Layout layout, typename idx_t>
idx_t geqrf_workquery(const Matrix<T, layout, idx_t>& A,
                      const Vector<T, idx_t>& tau);

/**
 * Compute the QR factorization of a general m-by-n matrix A.
 * This version of the function allocates the work array internally.
 */
template <typename T, Layout layout, typename idx_t, bool aligned = true>
void geqrf(const Matrix<T, layout, idx_t>& A,
           const Vector<T, idx_t>& tau)
{
    idx_t lwork = geqrf_workquery(A, tau);
    MemoryBlock<T, idx_t, aligned> work(lwork);
    geqrf(A, tau, work);
}

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_GEQRF_HPP