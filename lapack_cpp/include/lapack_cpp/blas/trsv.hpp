#ifndef LAPACK_CPP_TRSV_HPP
#define LAPACK_CPP_TRSV_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"
namespace lapack_cpp {

/**
 * Solve a triangular system of equations.
 */
template <typename T, Layout layout, typename idx_t>
void trsv(Uplo uplo,
          Op trans,
          Diag diag,
          const ConstMatrix<T, layout, idx_t>& A,
          const Vector<T, idx_t>& x);

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_TRSV_HPP