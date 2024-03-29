# LAPACKv4

This repository is a proposal for a new version of LAPACK: v4. The main change is that the main language of LAPACK would become C++ instead of Fortran.

LAPACK, written in Fortran in the 90s, and still heavily used today. It is difficult to understate the scope of LAPACK. Through years of development, it has grown in size to contain many different functionalities. There are the most popular factorizations, like the LU, the QR, SVD, and eigenvalue decompositions, but there are also endless variations of these routines.

Throughout the years, many different packages have popped up to smooth some of the rough edges of LAPACK. There are Fortran90 interfaces, C interfaces, C++ interfaces, Python interfaces, Julia interfaces, and even the Matlab platform started as an easy way to access LAPACK.

There are even libraries that go beyond just providing wrappers. LAPACKE, Eigen, Armadillo, Flame, SLATE, TLAPACK, (and probably many others I forget to mention) implement many of the common routines themselves and link to LAPACK for some of the more obscure routines, or sometimes they only implement basic versions of the routines, but provide the option to link to LAPACK if efficiency is required.

All this means that LAPACK still needs to be maintained, which is becoming increasingly difficult. Some problems include:

- Limited argument checking
- Limiting array access checking
- Multiple versions for different types: 4 types for real, double, complex, and complex double, and recently this got doubled to also include versions for 64 bit integers.
- Outdated testing framework
- Limited developer tools (Fortran is just not as popular as C/C++, it is much more difficult to find formatters, linters, and even programmers for Fortran)

The goal of this repo is to truly be a replacement for LAPACK, not just a wrapper on top of it. We believe a library that wishes to replace LAPACK should satisfy the following:

- Completeness, even the obscure or difficult routines should be implemented
- Full compatibility through Fortran API
- High performance (or at least as fast as the current LAPACK)

## Design principles

### 1. Make things easier for developers

The first, and most important design principle we have is that the library should make things easier for developers. This involves things like switching to C++, frequent assertions, and a good unit testing framework.

### 2. Matrix Storage

An important part of any linear algebra package is how they store matrices and vectors. In LAPACK, a matrix is stored in column-major order and passed by providing integers with its dimensions, a pointer to the data, and the leading dimension. While column-major order is the default in Fortran, row-major order is popular in C++. Many libraries decide to support both column- and row-major storage.

We have made the following decisions for matrices:
- Matrices are stored in a templated Matrix class. We believe that this is absolutely necessary to improve the developer experience. Since the Matrix class will know its size, it can easily check that no out-of-bounds accesses occur. This is even possible for submatrices that are locally defined, something that is not possible without such a class.
- Initially, we only support column-major storage. Our primary aim is to replace LAPACK, not provide a high-level interface that is easily callable (although a nice interface is preferable if easily achieved). Storing the matrices in column-major order makes it significantly easier to call LAPACK for any function that is yet to be translated. The Matrix classes do support row-major storage through an optional template argument so that we can easily support row-major storage at a later stage without an interface change.
- All matrices and vectors are non-owning lightweight wrappers. This allows easy wrapping of the pointers passed through the C/Fortran interface. We also need to use non-owning wrappers for submatrices anyway, so by making all matrices non-owning, we avoid having to support multiple types.
- An important part of a matrix class is constness. Typically in C++, if a non-owning wrapper is declared as const, the underlying data can still be modified. Our Matrix and Vector classes follow this rule and so somewhat counterintuitively, the following code does not complain: ```const Matrix<float> A(...); A(0,0) = 1.0;```. To declare a matrix as truly const, there is a separate class ```ConstMatrix``` which behaves the same as a normal matrix, except that the underlying data cannot be modified. This means that it can wrap a const pointer, only returns data by value, and returns ConstMatrices and ConstVectors for submatrices. A possible alternative would be to use ```Matrix<const float>```, however, this leads to the template parameter ```T``` being a const, which requires some nuance and we want to avoid issues for unaware programmers.


### 3. Interoperability with Fortran

We expect that it is infeasible to fully reimplement LAPACK in C++ or at least that this would require too much developer time to be worthwhile (remember that the entire point of the switch is to save development time). That is why we plan to design the library in a way that both Fortran and C++ can be mixed. For example, dgemv could have the following interfaces:

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
```c
void lapack_c_dgemv(char trans,
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

We keep all three interfaces intact at all times. Initially, this will be code written in Fortran, with C wrappers around it and C++ wrappers around the C wrappers. If a routine is rewritten in C++, we simply write C and Fortran wrappers around it. By developing this way, we can always keep the library in a "finished" state. Making it usable from the start. For relatively simple routines, like most of the BLAS routines, the wrappers will take more lines of code than the actual computational routines, but these are very simple wrappers that are easily implemented so their actual load on the programmer should not be too large.

Note: because we consider adding support for row-major storage at a later stage, it may be necessary to add a layout argument to the C interface initially so that supporting row-major storage does not require an interface change.

### 4. Future outlook

After translating all (or most of) the code to C++, it should be much easier to implement certain extra features. New, faster variants contributed by experts, support for row-major layout, ... 

### 5. Impact on current users

We plan to keep the impact on current users minimal, but there are still some changes to be expected:

- Installing lapackv4 will require compatible C++ and Fortran compilers, however, in the future, just a C++/C compiler could suffice.
- Our initial experiments have not shown a performance impact from mixing languages, but we cannot guarantee that this is the case for all compilers and use cases.

### 6. This repository

As a test for the interfaces, only ```gemm```, ```gemv```, and ```trsv``` have been translated.

```gemm``` is fully translated, and has a cpp implementation, with accompanying wrappers to call the code from C and Fortran. However, to allow the usage of optimized BLAS, it is possible to use disable this code and instead use C wrappers around Fortran, this is achieved using the flag ```USE_FORTRAN_BLAS```.

```gemv``` is not fully translated, only the wrappers to call the Fortran code from C and C++ are finished. Note that even if ```USE_FORTRAN_BLAS``` is false, ```gemv``` will call Fortran. The flag only has an effect on code that is translated to C++.

```trsv``` is also not fully translated. This tests an alternative way to define the c wrappers. Instead of messing around with fortran name mangling in c and then defining a separate wrapper, the c wrapper is defined in fortran. This requires less code and is the "proper" way to mix languages according to the gfortran documentation.

Some examples can be found in the ```example/``` folder.