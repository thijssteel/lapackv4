#ifndef LAPACK_CPP_TRSV_HPP
#define LAPACK_CPP_TRSV_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"
namespace lapack_cpp {

/**
 * Solve a triangular system of equations.
 */
template <typename T, typename idx_t, Layout layout>
void trsv(Uplo uplo,
          Op trans,
          Diag diag,
          const ConstMatrix<T, idx_t, layout>& A,
          const Vector<T, idx_t>& x);

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_TRSV_HPP