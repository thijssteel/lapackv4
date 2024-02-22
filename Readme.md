# LAPACKv4

This repository is a proposal for a new version of LAPACK: v4. The main change is that the main language of LAPACK would become C++ instead of Fortran.

LAPACK, written in Fortran in the 90s and still heavily used today. It difficult to understate the scope of LAPACK. Through years of development, it has grown in size to contain many different functionalities. There are the most popular factorizations, like the LU, the QR, SVD and eigenvalue decompositions, but there are also endless variations of these routines.

Throughout the years, many different packages have popped up to smooth some of the rough edges of LAPACK. There are Fortran90 interfaces, C interfaces, C++ interfaces, Python interfaces, Julia interfaces and even the Matlab platform started as an easy way to access LAPACK.

There are even libraries which go beyond just providing wrappers. LAPACKE, Eigen, Armadillo, Flame, SLATE, TLAPACK, (and probably many others I forget to mention) implement many of the common routines themselves and link to LAPACK for some of the more obscure routines, or sometimes they only implement basic versions of the routines, but provide the option to link to LAPACK if efficiency is required.

All this means that LAPACK still needs to be maintained, which is becoming increasingly difficult. Some problems include:

- Limited argument checking
- Limiting array access checking
- Multiple versions for different types: 4 types for real, double, complex and complex double, and recently this got double to also include versions for 64 bit integers.
- Outdated testing framework
- Limited developer tools (Fortran is just not as popular as C/C++, it is much more difficult to find formatters, linters, and even programmers for Fortran)

The goal of this repo is to truly be a replacement for LAPACK, not just a wrapper on top of it. We believe a library that wishes to replace LAPACK should satisfy the following:

- Completeness, even the obscure or difficult routines should be implemented
- Full compatibility through Fortran api
- High performance (or at least as fast as the current LAPACK)

## Design principles

### 1. Make things easier for developers

The first, and most important design principle we have is that the library should make things easier for developers. This involves things like switching to C++, frequent assertions, and a good unit testing framework.

### 2. Matrix storage

An important part of any linear algebra package is how they store matrices and vectors. In LAPACK, matrices are stored in column-major order and passed by providing integers with its dimensions, a pointer to the data and the leading dimension. While column-major order is the default in Fortran, row-major order is popular in C++. Many libraries decide to support both column- and row-major storage.

We have made the following decisions for matrices:
- Matrices are stored in a templated Matrix class. We believe that this is absolutely necessary to improve the developer experience. Since the Matrix class will know its size, it can easy check that no out of bounds accesses occur. This is even possible for submatrices that are locally defined, something that is not possible without such a class.
- Initially, we only support column-major order. Our primary aim is to replace LAPACK, not provide a high-level interface that is easily callable (although a nice interface is preferable if easily achieved). Storing the matrices in column-major order makes it significantly easier to call LAPACK for any function that is yet to be translated. The Matrix classes do support row-major storage through an optional template argument.
- All matrices and vectors are non-owning lightweight wrappers. This allows easy wrapping of the pointers passed through the C/Fortran interface. We also need to use non-owning wrappers for submatrices anyway, so by making all matrices non-owning, we avoid having to support multiple types.
- An important part of a matrix class is constness. Typically in C++, if a non-owning wrapper is declared as const, the underlying data can still be modified. Our Matrix and Vector classes follow this rule and so somewhat counterintuitively, the following code does not complain: ```const Matrix<float> A(...); A(0,0) = 1.0;```. To declare a matrix as truly const, there is a separate class ```ConstMatrix``` which behaves the same as a normal matrix, except that the underlying data cannot be modified. This means that it can wrap a const pointer, only returns data by value and return ConstMatrices and ConstVectors for submatrices. A possible alternative would be to use ```Matrix<const float>```, however, this leads to the template parameter ```T``` being a const, which requires some nuance and we want to avoid issues for unaware programmers.


### 3. Interoperability with Fortran

We expect that it is infeasible to fully reimplement LAPACK in C++ or at least that this would require too much developer time to be worthwhile (remember that the entire point of the switch is to save development time). That is why we plan to design the library in a way that both fortran and c++ can be mixed in a consistent way. For example, dgemv could have the following interfaces:

In C++:
```c++
template <typename T, typename idx_t>
void gemv(Op trans,
        T alpha,
        const ConstMatrix<T, idx_t>& A,
        const ConstVector<T, idx_t>& x,
        T beta,
        const Vector<T, idx_t>& y)
```
In C:
```c++
void dgemv_c(char trans,
                int m,
                int n,
                double alpha,
                const double* A,
                int lda,
                const double* x,
                int incx,
                double beta,
                double* y,
                int incy)
```
In Fortran it would have the traditional interface:
```Fortran
subroutine dgemv( trans, m, n, alpha, A, lda, x, incx, beta, y, incy )
```

We keep all three interfaces intact at all times. Initially, this will be code written in Fortran, with C wrappers around it and C++ wrappers around those. If a routine is rewritten in C++, we simply write C and Fortran wrappers around it. By developing this way, we can always keep the library in a "finished" state. Making it usable from the start.

Note: because we consider adding support for row-major layout at a later stage, it may be necessary to add a layout argument to the C interface initially so that supporting row-major layout does not require an interface change.

### 4. Future outlook

After translating all (or most of) the code to C++, it should be much easier to implement certain extra features. New, faster variants contributed by experts, support for row-major layout, ... 

### 5. Impact on current users

We plan to keep the impact on current users minimal, but there are still some changes to be expected:

- Installing lapackv4 will require compatible C++ and Fortran compilers, however, in the future, just a C++/C compiler could suffice.
- Our initial experiments have not shown a performance impact from mixing languages, but we cannot guarantee that this is the case for all compilers and use cases.

### 6. This repository

As a test for the interfaces, only the general matrix-matrix multiplication has been translated (gemm).

There are two examples: `test_gemm.cpp` and `test_gemm.f90`. The provided makefile will compile these into 4 files
- test_gemm_fortran_fortran: example calls pure fortran
- test_gemm_fortran_cpp: example calls cpp code through fortran interface
- test_gemm_cpp_cpp: example calls pure cpp
- test_gemm_cpp_fortran: example calls fortran code through cpp interface. Note, it is linked with BLAS via -lblas instead of the code in this repo. This is to illustrate the purpose of this interface. It can be used initially while there is no C++ implementation available or later when we want to link with optimized blas. Vendors could even provide their own cpp implementation.

