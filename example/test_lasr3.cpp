
#include <complex>
#include <iostream>

#include "lapack_cpp.hpp"

using namespace lapack_cpp;

template <typename T,
          Layout layout = Layout::ColMajor,
          typename idx_t = lapack_idx_t>
void test_lasr3(idx_t m, idx_t n, idx_t k, Side side, Direction direction)
{
    idx_t n_rot = side == Side::Left ? m - 1 : n - 1;

    MemoryBlock<T, idx_t> A_(m, n);
    Matrix<T, layout, idx_t> A(m, n, A_);

    MemoryBlock<real_t<T>, idx_t> C_(n_rot, k);
    Matrix<real_t<T>, layout, idx_t> C(n_rot, k, C_);

    MemoryBlock<T, idx_t> S_(n_rot, k);
    Matrix<T, layout, idx_t> S(n_rot, k, S_);

    randomize(A);

    // Generate random, but valid rotations
    randomize(C);
    randomize(S);
    for (idx_t i = 0; i < n_rot; ++i) {
        for (idx_t j = 0; j < k; ++j) {
            T f = C(i, j);
            T g = S(i, j);
            T r;
            lartg(f, g, C(i, j), S(i, j), r);
        }
    }

    MemoryBlock<T, idx_t> A_copy_(m, n);
    Matrix<T, layout, idx_t> A_copy(m, n, A_copy_);

    A_copy = A;

    // Apply rotations using simple loops as test
    if (side == Side::Left) {
        if (direction == Direction::Forward) {
            for (idx_t i = 0; i < k; ++i)
                for (idx_t j = 0; j < n_rot; ++j)
                    rot(A_copy.row(j), A_copy.row(j + 1), C(j, i), S(j, i));
        }
        else {
            for (idx_t i = 0; i < k; ++i)
                for (idx_t j = n_rot - 1; j >= 0; --j)
                    rot(A_copy.row(j), A_copy.row(j + 1), C(j, i), S(j, i));
        }
    }
    else {
        if (direction == Direction::Forward) {
            for (idx_t i = 0; i < k; ++i)
                for (idx_t j = 0; j < n_rot; ++j)
                    rot(A_copy.column(j), A_copy.column(j + 1), C(j, i),
                        S(j, i));
        }
        else {
            for (idx_t i = 0; i < k; ++i)
                for (idx_t j = n_rot - 1; j >= 0; --j)
                    rot(A_copy.column(j), A_copy.column(j + 1), C(j, i),
                        S(j, i));
        }
    }

    MemoryBlock<T, idx_t, true> work(k + 2, n, Layout::RowMajor);

    lasr3(side, direction, C.as_const(), S.as_const(), A, work);

    // Check that the result is the same
    for (idx_t i = 0; i < m; ++i)
        for (idx_t j = 0; j < n; ++j)
            A_copy(i, j) -= A(i, j);

    print(A_copy);
}

int main()
{
    typedef lapack_idx_t idx_t;

    for (idx_t n_rot = 3; n_rot < 10; ++n_rot) {
        for (idx_t k = 1; k <= 2; k++) {
            std::cout << "n_rot = " << n_rot << ", k = " << k << std::endl;
            test_lasr3<float, Layout::ColMajor, lapack_idx_t>(
                n_rot + 1, 1, k, Side::Left, Direction::Forward);
            test_lasr3<float, Layout::ColMajor, lapack_idx_t>(
                1, n_rot + 1, k, Side::Right, Direction::Forward);
            test_lasr3<float, Layout::ColMajor, lapack_idx_t>(
                n_rot + 1, 1, k, Side::Left, Direction::Backward);
            test_lasr3<float, Layout::ColMajor, lapack_idx_t>(
                1, n_rot + 1, k, Side::Right, Direction::Backward);
        }
    }
    return 0;
}