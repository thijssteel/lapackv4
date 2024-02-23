#ifndef LAPACK_CPP_ROT_HPP
#define LAPACK_CPP_ROT_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"
namespace lapack_cpp {

/**
 * Apply a plane rotation to a pair of vectors.
 */
template <typename T, typename TC, typename TS, typename idx_t>
void rot(const Vector<T, idx_t>& x, const Vector<T, idx_t>& y, TC c, TS s);

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_GEMV_HPP