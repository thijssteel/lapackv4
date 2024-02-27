#ifndef LAPACK_CPP_BASE_TYPES_HPP
#define LAPACK_CPP_BASE_TYPES_HPP

#include <cstdint>
#include <complex>
#include <type_traits>

// While the library is templated and can handle different integers,
// this is the default integer type used by the library.
typedef int64_t lapack_idx_t;

// declare conj for real types
template <typename T,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
inline T conj(const T& x)
{
    return x;
}

// declare real for real types
template <typename T,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
inline T real(const T& x)
{
    return x;
}


#endif // LAPACK_CPP_BASE_TYPES_HPP