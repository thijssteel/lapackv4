CXX=g++

CXXFLAGS=-I$(incdir) -std=c++20 -Wall -Wextra -pedantic -g

all: ./example/test_gemm_fortran_fortran ./example/test_gemm_fortran_cpp ./example/test_gemm_cpp_cpp ./example/test_gemm_cpp_fortran

# Fortran files

lapack_fortran/src/%.o: lapack_fortran/src/%.f
	gfortran -c -o $@ $<

lapack_fortran/cwrappers/%.o: lapack_fortran/cwrappers/%.f90
	gfortran -c -o $@ $<

lapack_fortran/cwrappers/dgemm.o: lapack_fortran/cwrappers/kinds.o
lapack_fortran/cwrappers/sgemm.o: lapack_fortran/cwrappers/kinds.o

# C files

lapack_c/cppwrappers/%.o: lapack_c/cppwrappers/%.cpp
	g++ -std=c++17 -I./lapack_c/include -I./lapack_cpp/include -c -o $@ $<

lapack_c/fortranwrappers/%.o: lapack_c/fortranwrappers/%.c
	gcc -I./lapack_c/include -I./lapack_cpp/include -c -o $@ $<

# C++ files
lapack_cpp/cwrappers/blas/%.o: lapack_cpp/cwrappers/blas/%.cpp
	g++ -std=c++17 -I./lapack_c/include -I./lapack_cpp/include  -c -o $@ $<

# Test files
# The compilation commands here are very important, depending on what files are being linked to, we use the C++ or Fortran implementation.
# Also important, if only C++ implementations are being used, we don't need a fortran compiler.

example/test_gemm_fortran_fortran: example/test_gemm.f90 lapack_fortran/src/sgemm.o lapack_fortran/src/lsame.o  lapack_fortran/src/xerbla.o
	gfortran -o $@ $^

example/test_gemm_fortran_cpp: example/test_gemm.f90 lapack_fortran/cwrappers/sgemm.o lapack_c/cppwrappers/gemm.o
	gfortran -o $@ $^ -lstdc++

example/test_gemm_cpp_cpp: example/test_gemm.cpp
	g++ -std=c++17 -I./lapack_cpp/include -o $@ $<

example/test_gemm_cpp_fortran: example/test_gemm.cpp lapack_cpp/cwrappers/blas/gemm.o lapack_c/fortranwrappers/gemm.o
	g++ -std=c++17 -I./lapack_c/include/ -I./lapack_cpp/include/ -o $@ $^ -lgfortran -lblas

clean:
	find . -type f -name '*.o' -delete