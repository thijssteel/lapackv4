#ifndef LAPACK_CPP_GEMV_HPP
#define LAPACK_CPP_GEMV_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"
namespace lapack_cpp {

/**
 * Perform a matrix-vector multiplication.
 */
template <typename T, Layout layout, typename idx_t>
void gemv(Op trans,
          T alpha,
          const ConstMatrix<T, layout, idx_t>& A,
          const ConstVector<T, idx_t>& x,
          T beta,
          const Vector<T, idx_t>& y);

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_GEMV_HPP