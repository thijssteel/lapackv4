
#include <iostream>
#include <complex>

#include "lapack_cpp.hpp"

using namespace lapack_cpp;

template <typename T>
void test_blas(){
    int m = 5;
    int n = 3;

    MemoryBlock<T> A_(m,n);
    Matrix<T> A(m,n, A_);

    MemoryBlock<T> tau_(n);
    Vector<T> tau(n, tau_);

    randomize(A);

    std::cout << "A before QR = " << std::endl;
    print(A);
    
    geqrf(A, tau);

    std::cout << "A after QR = " << std::endl;
    print(A);
}



int main(){

    test_blas<float>();
    test_blas<std::complex<float>>();
    test_blas<double>();
    test_blas<std::complex<double>>();

    return 0;
}