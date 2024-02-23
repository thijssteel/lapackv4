
#include <iostream>
#include <complex>

#include "lapack_cpp.hpp"

using namespace lapack_cpp;

template <typename T>
void test_blas(){
    int m = 2;
    int n = 3;

    MemoryBlock<T> A_(m,n);
    Matrix<T> A(m,n, A_);

    MemoryBlock<T> x_(n);
    Vector<T> x(n, x_);

    MemoryBlock<T> y_(m);
    Vector<T> y(m, y_);

    randomize(A);
    randomize(x);

    // Calculate y = A*x with loops
    for (int i = 0; i < m; ++i){
        y[i] = 0;
        for (int j = 0; j < n; ++j){
            y[i] += A(i,j)*x[j];
        }
    }

    gemv(Op::NoTrans, (T) 1.0, A.as_const(), x.as_const(), (T) -1., y);

    print(y);
}



int main(){


    test_blas<float>();
    test_blas<std::complex<float>>();
    test_blas<double>();
    test_blas<std::complex<double>>();


    return 0;
}