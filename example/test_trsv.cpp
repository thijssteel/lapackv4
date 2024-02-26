
#include <iostream>
#include <complex>

#include "lapack_cpp.hpp"

using namespace lapack_cpp;

template <typename T>
void test_blas(){
    int n = 5;

    MemoryBlock<T> A_(n,n);
    Matrix<T> A(n,n, A_);

    MemoryBlock<T> x_(n);
    Vector<T> x(n, x_);

    randomize(A);

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            if (i < j){
                A(i,j) = 0;
            }
        }
    }

    print(A);

    randomize(x);

    print(x);

    lapack_cpp::trsv(Uplo::Lower, Op::NoTrans, Diag::NonUnit, A.as_const(), x);

    print(x);
}



int main(){


    test_blas<float>();
    test_blas<std::complex<float>>();
    test_blas<double>();
    test_blas<std::complex<double>>();


    return 0;
}