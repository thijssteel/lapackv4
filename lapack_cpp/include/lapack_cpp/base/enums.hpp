#ifndef LAPACK_CPP_ENUMP_HPP
#define LAPACK_CPP_ENUMP_HPP

#include <cassert>

namespace lapack_cpp {

enum class Op : char { NoTrans = 'N', Trans = 'T', ConjTrans = 'C' };

inline constexpr Op char2op(char t)
{
    switch (t) {
        case 'N':
            return Op::NoTrans;
        case 'T':
            return Op::Trans;
        case 'C':
            return Op::ConjTrans;
        default:
            assert(false);
    }
}

enum class Side : char { Left = 'L', Right = 'R' };

inline constexpr Side char2side(char t)
{
    switch (t) {
        case 'L':
            return Side::Left;
        case 'R':
            return Side::Right;
        default:
            assert(false);
    }
}

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_ENUMP_HPP