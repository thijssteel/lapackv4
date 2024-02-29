#ifndef LAPACK_CPP_LASR3_HPP
#define LAPACK_CPP_LASR3_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/blas/rot.hpp"
#include "lapack_cpp/utils.hpp"

namespace lapack_cpp {

/**
 * Applies a sequence of plane rotations to a general rectangular matrix.
 *
 * @param side   Specifies whether the plane rotations are applied to A on the
 *               left or right.
 *
 * @param direct Specifies whether the sequence is applied in forward or
 *               backward order.
 *
 * @param C      The cosines of the rotations
 *               If side == Side::Left, C is an array of dimension m-1 x k
 *               If side == Side::Right, C is an array of dimension n-1 x k
 *
 * @param S      The sines of the rotations
 *               If side == Side::Left, S is an array of dimension m-1 x k
 *               If side == Side::Right, S is an array of dimension n-1 x k
 *
 * @param A      m x n matrix to which the rotations are to be applied.
 *
 * @param work   A workspace array, its dimension is specified by lasr3_work.
 *               Note, because we use the work array to align the data for
 *               efficient vectorization, this routine only support aligned
 *               work arrays.
 *
 * @tparam T     The type of the elements of the matrix A.
 *
 * @tparam TC    The type of the elements of the matrix C.
 *
 * @tparam TS    The type of the elements of the matrix S.
 *
 * @tparam layout The layout of the matrices A, C, and S.
 *
 * @tparam idx_t  The type of the indices.
 *
 * @note if either TC or TS is a complex type, then T must also be a complex
 * type.
 *
 */
template <typename T, typename TC, typename TS, Layout layout, typename idx_t>
void lasr3(Side side,
           Direction direct,
           ConstMatrix<TC, layout, idx_t> C,
           ConstMatrix<TS, layout, idx_t> S,
           Matrix<T, layout, idx_t> A,
           MemoryBlock<T, idx_t, true> work);

/** Applies two plane rotations to three vectors.
 *
 * [x1]   [1 0   0  ]   [c1  s1 0]   [x1]
 * [x2] = [0 c2  s2 ] * [-s1 c1 0] * [x2]
 * [x3]   [0 -s2 c2 ]   [0   0  1]   [x3]
 *
 * @param x1  The first vector.
 * @param x2  The second vector.
 * @param x3  The third vector.
 * @param c1  The cosine of the first rotation.
 * @param s1  The sine of the first rotation.
 * @param c2  The cosine of the second rotation.
 * @param s2  The sine of the second rotation.
 */
template <typename T, typename TC, typename TS, typename idx_t>
void rot_fuse2x1(Vector<T, idx_t> x1,
                 Vector<T, idx_t> x2,
                 Vector<T, idx_t> x3,
                 TC c1,
                 TS s1,
                 TC c2,
                 TS s2)
{
    assert(x1.size() == x2.size() && x2.size() == x3.size());
    const idx_t n = x1.size();
    for (idx_t i = 0; i < n; ++i) {
        T x2_g1 = -conj(s1) * x1[i] + c1 * x2[i];
        x1[i] = c1 * x1[i] + s1 * x2[i];
        x2[i] = c2 * x2_g1 + s2 * x3[i];
        x3[i] = -conj(s2) * x2_g1 + c2 * x3[i];
    }
}

/** Applies two plane rotations to three vectors.
 *
 * [x1]   [c2  s2 0]   [1 0   0  ]   [x1]
 * [x2] = [-s2 c2 0] * [0 c1  s1 ] * [x2]
 * [x3]   [0   0  1]   [0 -s1 c1 ]   [x3]
 *
 * @param x1  The first vector.
 * @param x2  The second vector.
 * @param x3  The third vector.
 * @param c1  The cosine of the first rotation.
 * @param s1  The sine of the first rotation.
 * @param c2  The cosine of the second rotation.
 * @param s2  The sine of the second rotation.
 */
template <typename T, typename TC, typename TS, typename idx_t>
void rot_fuse1x2(Vector<T, idx_t> x1,
                 Vector<T, idx_t> x2,
                 Vector<T, idx_t> x3,
                 TC c1,
                 TS s1,
                 TC c2,
                 TS s2)
{
    assert(x1.size() == x2.size() && x2.size() == x3.size());
    const idx_t n = x1.size();
    for (idx_t i = 0; i < n; ++i) {
        T x2_g1 = c1 * x2[i] + s1 * x3[i];
        x3[i] = -conj(s1) * x2[i] + c1 * x3[i];
        x2[i] = -conj(s2) * x1[i] + c2 * x2_g1;
        x1[i] = c2 * x1[i] + s2 * x2_g1;
    }
}

/** Applies four plane rotations to four vectors.
 *
 * [x1]   [1 0   0  0]   [1 0 0   0 ]   [c2  s2 0 0]   [1 0   0  0]   [x1]
 * [x2] = [0 c4  s4 0] * [0 1 0   0 ] * [-s2 c2 0 0] * [0 c1  s1 0] * [x2]
 * [x3]   [0 -s4 c4 0]   [0 0 c3  s3]   [0   0  1 0]   [0 -s1 c1 0]   [x3]
 * [x4]   [0 0   0  1]   [0 0 -s3 c3]   [0   0  0 1]   [0 0   0  0]   [x4]
 *
 * Note: the rotations are applied in the order G1, G2, G3, G4,
 * but the order G1, G3, G2, G4 is also possible.
 *
 * @param x1  The first vector.
 * @param x2  The second vector.
 * @param x3  The third vector.
 * @param x4  The fourth vector.
 * @param c1  The cosine of the first rotation.
 * @param s1  The sine of the first rotation.
 * @param c2  The cosine of the second rotation.
 * @param s2  The sine of the second rotation.
 * @param c3  The cosine of the third rotation.
 * @param s3  The sine of the third rotation.
 * @param c4  The cosine of the fourth rotation.
 * @param s4  The sine of the fourth rotation.
 */
template <typename T, typename TC, typename TS, typename idx_t>
void rot_fuse2x2(Vector<T, idx_t> x1,
                 Vector<T, idx_t> x2,
                 Vector<T, idx_t> x3,
                 Vector<T, idx_t> x4,
                 TC c1,
                 TS s1,
                 TC c2,
                 TS s2,
                 TC c3,
                 TS s3,
                 TC c4,
                 TS s4)
{
    assert(x1.size() == x2.size() && x1.size() == x3.size() &&
           x1.size() == x4.size());
    const idx_t n = x1.size();
    for (idx_t i = 0; i < n; ++i) {
        T x2_g1 = c1 * x2[i] + s1 * x3[i];
        T x3_g1 = -conj(s1) * x2[i] + c1 * x3[i];
        T x2_g2 = -conj(s2) * x1[i] + c2 * x2_g1;
        x1[i] = c2 * x1[i] + s2 * x2_g1;
        T x3_g3 = c3 * x3_g1 + s3 * x4[i];
        x4[i] = -conj(s3) * x3_g1 + c3 * x4[i];
        x2[i] = c4 * x2_g2 + s4 * x3_g3;
        x3[i] = -conj(s4) * x2_g2 + c4 * x3_g3;
    }
}

// Kernel for lasr3, forward left variant
template <typename T, typename TC, typename TS, Layout layout, typename idx_t>
void lasr3_kernel_forward_left(ConstMatrix<TC, layout, idx_t> C,
                               ConstMatrix<TS, layout, idx_t> S,
                               Matrix<T, layout, idx_t> A,
                               MemoryBlock<T, idx_t, true> work)
{
    idx_t m = A.num_rows();
    idx_t n = A.num_columns();
    idx_t k = C.num_columns();

    assert(C.num_rows() == S.num_rows());
    assert(C.num_columns() == S.num_columns());
    assert(m == C.num_rows() + 1);
    assert(m > k + 1);

    // Number of rows that will be stored in the packed workspace matrix
    const idx_t np = k + 2;
    // During the algorithm, instead of applying the rotations directly to
    // the matrix A, we will apply them to a packed version of A. This is
    // done to allow for efficient vectorization of the rotations. The
    // packed matrix is stored in the workspace array. Since we apply the
    // rotations from the left, this packed matrix is always row-major
    // regardless of the layout of A.
    Matrix<T, Layout::RowMajor, idx_t> A_pack(np, n, work);
    // A_pack will store np rows of A.
    // To avoid having to shuffle rows in A, we keep track of two indices
    // For i = i_pack2, ..., i_pack2 + np - 1, row (i_pack + i - i_pack2 +
    // np) % np of A_pack corresponds to row i of A2.
    idx_t i_pack = 0;
    idx_t i_pack2 = 0;

    // Copy the first np rows of A to A_pack
    for (idx_t i = 0; i < np; ++i) {
        for (idx_t j = 0; j < n; ++j) {
            A_pack(i, j) = A(i, j);
        }
    }

    // Startup phase
    for (idx_t j = 0; j + 1 < k; ++j) {
        for (idx_t i = 0, g = j; i < j + 1; ++i, --g) {
            rot(A_pack.row(g), A_pack.row(g + 1), C(g, i), S(g, i));
        }
    }

    // Pipeline phase
    for (idx_t j = k - 1; j + 1 < m - 1; j += 2) {
        for (idx_t i = 0, g = j; i + 1 < k; i += 2, g -= 2) {
            auto a1 = A_pack.row((g - 1 + i_pack - i_pack2 + np) % np);
            auto a2 = A_pack.row((g + i_pack - i_pack2 + np) % np);
            auto a3 = A_pack.row((g + 1 + i_pack - i_pack2 + np) % np);
            auto a4 = A_pack.row((g + 2 + i_pack - i_pack2 + np) % np);

            rot_fuse2x2(a1, a2, a3, a4, C(g, i), S(g, i), C(g - 1, i + 1),
                        S(g - 1, i + 1), C(g + 1, i), S(g + 1, i), C(g, i + 1),
                        S(g, i + 1));
        }
        if (k % 2 == 1) {
            // k is odd, so there are two more rotations to apply
            idx_t i = k - 1;
            idx_t g = j - i;

            auto a1 = A_pack.row((g + i_pack - i_pack2 + np) % np);
            auto a2 = A_pack.row((g + 1 + i_pack - i_pack2 + np) % np);
            auto a3 = A_pack.row((g + 2 + i_pack - i_pack2 + np) % np);

            rot_fuse2x1(a1, a2, a3, C(g, i), S(g, i), C(g + 1, i), S(g + 1, i));
        }
        // rows i_pack and i_pack+1 of A_pack are finished, copy them back
        // to A
        for (idx_t i = 0; i < n; i++) {
            A(i_pack2, i) = A_pack(i_pack, i);
            A(i_pack2 + 1, i) = A_pack((i_pack + 1) % np, i);
        }
        // Pack next rows and update i_pack and i_pack2
        if (j + 4 < m) {
            for (idx_t i = 0; i < n; i++) {
                A_pack(i_pack, i) = A(j + 3, i);
                A_pack((i_pack + 1) % np, i) = A(j + 4, i);
            }
            i_pack = (i_pack + 2) % np;
            i_pack2 += 2;
        }
        else if (j + 3 < m) {
            for (idx_t i = 0; i < n; i++) {
                A_pack(i_pack, i) = A(j + 3, i);
            }
            i_pack = (i_pack + 1) % np;
            i_pack2 += 1;
        }
    }

    // Shutdown phase
    for (idx_t j = (m - k + 1) % 2; j < k; ++j) {
        for (idx_t i = j, g = m - 2; i < k; ++i, --g) {
            auto a1 = A_pack.row((g + i_pack - i_pack2 + np) % np);
            auto a2 = A_pack.row((g + 1 + i_pack - i_pack2 + np) % np);
            rot(a1, a2, C(g, i), S(g, i));
        }
    }

    // Copy the last np rows of A_pack back to A
    for (idx_t i = 0; i < std::min(np, m - i_pack2); ++i) {
        for (idx_t j = 0; j < n; ++j) {
            A(i_pack2 + i, j) = A_pack((i_pack + i) % np, j);
        }
    }
}

// Kernel for lasr3, forward left variant
template <typename T, typename TC, typename TS, Layout layout, typename idx_t>
void lasr3_kernel_backward_left(ConstMatrix<TC, layout, idx_t> C,
                                ConstMatrix<TS, layout, idx_t> S,
                                Matrix<T, layout, idx_t> A,
                                MemoryBlock<T, idx_t, true> work)
{
    idx_t m = A.num_rows();
    idx_t n = A.num_columns();
    idx_t k = C.num_columns();

    assert(C.num_rows() == S.num_rows());
    assert(C.num_columns() == S.num_columns());
    assert(m == C.num_rows() + 1);
    assert(m > k + 1);

    // Number of rows that will be stored in the packed workspace matrix
    const idx_t np = k + 2;
    // During the algorithm, instead of applying the rotations directly to
    // the matrix A, we will apply them to a packed version of A. This is
    // done to allow for efficient vectorization of the rotations. The
    // packed matrix is stored in the workspace array. Since we apply the
    // rotations from the left, this packed matrix is always row-major
    // regardless of the layout of A.
    Matrix<T, Layout::RowMajor, idx_t> A_pack(np, n, work);
    // A_pack will store np rows of A.
    // To avoid having to shuffle rows in A, we keep track of two indices
    // For i = i_pack2, ..., i_pack2 + np - 1, row (i_pack + i - i_pack2 +
    // np) % np of A_pack corresponds to row i of A2.
    idx_t i_pack = 0;
    idx_t i_pack2 = m - np;

    // Copy the last np rows of A to A_pack
    for (idx_t i = 0; i < np; ++i) {
        for (idx_t j = 0; j < n; ++j) {
            A_pack(i, j) = A(m - np + i, j);
        }
    }

    // Startup phase
    for (idx_t j = 0; j + 1 < k; ++j) {
        for (idx_t i = 0, g = m - 2 - j; i < j + 1; ++i, ++g) {
            rot(A_pack.row((g - i_pack2) % np),
                A_pack.row((g + 1 - i_pack2) % np), C(g, i), S(g, i));
        }
    }

    // Pipeline phase
    for (idx_t j = k - 1; j + 1 < m - 1; j += 2) {
        for (idx_t i = 0, g = m - 2 - j; i + 1 < k; i += 2, g += 2) {
            auto a1 = A_pack.row((g - 1 + i_pack - i_pack2) % np);
            auto a2 = A_pack.row((g + i_pack - i_pack2) % np);
            auto a3 = A_pack.row((g + 1 + i_pack - i_pack2) % np);
            auto a4 = A_pack.row((g + 2 + i_pack - i_pack2) % np);

            rot_fuse2x2(a1, a2, a3, a4, C(g, i), S(g, i), C(g - 1, i),
                        S(g - 1, i), C(g + 1, i + 1), S(g + 1, i + 1),
                        C(g, i + 1), S(g, i + 1));
        }
        if (k % 2 == 1) {
            // k is odd, so there are two more rotations to apply
            idx_t i = k - 1;
            idx_t g = m - 2 - j + i;

            auto a1 = A_pack.row((g - 1 + i_pack - i_pack2 + np) % np);
            auto a2 = A_pack.row((g + i_pack - i_pack2 + np) % np);
            auto a3 = A_pack.row((g + 1 + i_pack - i_pack2 + np) % np);

            rot_fuse1x2(a1, a2, a3, C(g, i), S(g, i), C(g - 1, i), S(g - 1, i));
        }
        // rows i_pack+np-2 and i_pack+np-1 are finished, copy them back to A
        for (idx_t i = 0; i < n; i++) {
            A(i_pack2 + np - 2, i) = A_pack((i_pack + np - 2) % np, i);
            A(i_pack2 + np - 1, i) = A_pack((i_pack + np - 1) % np, i);
        }
        // Pack next rows and update i_pack and i_pack2
        if (i_pack2 > 1) {
            for (idx_t i = 0; i < n; i++) {
                A_pack((i_pack + np - 2) % np, i) = A(i_pack2 - 2, i);
                A_pack((i_pack + np - 1) % np, i) = A(i_pack2 - 1, i);
            }
            i_pack = (i_pack + np - 2) % np;
            i_pack2 -= 2;
        }
        else if (i_pack2 > 0) {
            for (idx_t i = 0; i < n; i++) {
                A_pack((i_pack + np - 1) % np, i) = A(i_pack2 - 1, i);
            }
            i_pack = (i_pack + np - 1) % np;
            i_pack2 -= 1;
        }
    }

    // By now, i_pack2 should point to the start of the matrix.
    assert(i_pack2 == 0);

    // Shutdown phase
    for (idx_t j = (m - k + 1) % 2; j < k; ++j) {
        for (idx_t i = j, g = 0; i < k; ++i, ++g) {
            auto a1 = A_pack.row((g + i_pack - i_pack2 + np) % np);
            auto a2 = A_pack.row((g + 1 + i_pack - i_pack2 + np) % np);
            rot(a1, a2, C(g, i), S(g, i));
        }
    }

    // Copy the first np rows of A_pack back to A
    for (idx_t i = 0; i < np; ++i) {
        for (idx_t j = 0; j < n; ++j) {
            A(i, j) = A_pack((i_pack + i) % np, j);
        }
    }
}

/**
 * @copydoc lasr3
 */
template <typename T, typename TC, typename TS, Layout layout, typename idx_t>
void lasr3(Side side,
           Direction direct,
           ConstMatrix<TC, layout, idx_t> C,
           ConstMatrix<TS, layout, idx_t> S,
           Matrix<T, layout, idx_t> A,
           MemoryBlock<T, idx_t, true> work)
{
    idx_t m = A.num_rows();
    idx_t n = A.num_columns();
    idx_t k = C.num_columns();

    assert(side == Side::Left || side == Side::Right);
    assert(direct == Direction::Forward || direct == Direction::Backward);
    assert(C.num_rows() == S.num_rows());
    assert(C.num_columns() == S.num_columns());
    assert(side == Side::Left ? m : n == C.num_rows() + 1);
    assert((side == Side::Left ? m : n) > k + 1);

    if (side == Side::Left and direct == Direction::Forward) {
        // TODO: do blocking
        lasr3_kernel_forward_left(C, S, A, work);
    }

    if (side == Side::Left and direct == Direction::Backward) {
        // TODO: do blocking
        lasr3_kernel_backward_left(C, S, A, work);
    }
}

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_LASR3_HPP