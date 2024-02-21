
#include <iostream>
#include <complex>

#include "lapack_cpp.hpp"


using namespace lapack_cpp;

template <typename T>
void test_blas(){
    int m = 3;
    int n = 2;
    int k = 1;

    MemoryBlock<T> A_(m,k);
    Matrix<T> A(m,k, A_);

    MemoryBlock<T> B_(k,n);
    Matrix<T> B(k, n, B_);

    MemoryBlock<T> C_(m,n);
    Matrix<T> C(m, n, C_);

    randomize(A);
    randomize(B);

    // Calculate C = A*B with loops
    for (int i = 0; i < m; ++i){
        for (int j = 0; j < n; ++j){
            C(i,j) = 0;
            for (int l = 0; l < k; ++l){
                C(i,j) += A(i,l)*B(l,j);
            }
        }
    }

    gemm(Op::NoTrans, Op::NoTrans, (T) 1.0, A.as_const(), B.as_const(), (T) -1., C);

    print(C);
}



int main(){


    test_blas<float>();
    test_blas<std::complex<float>>();
    test_blas<double>();
    test_blas<std::complex<double>>();


    return 0;
}