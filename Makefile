CXX=g++

CXXFLAGS=-I$(incdir) -std=c++20 -Wall -Wextra -pedantic -g

all: ./bin/test_matrix ./bin/test_gemm ./bin/test_gemv ./bin/test_gemv_fort ./bin/test_geqrf

# Fortran files

lapack_fortran/src/*.o: lapack_fortran/src/*.f
	gfortran -c -o $@ $<

lapack_fortran/cwrappers/kinds.o: lapack_fortran/cwrappers/kinds.f90
	gfortran -c -o $@ $<

lapack_fortran/cwrappers/sgemm.o: lapack_fortran/cwrappers/sgemm.f90
	gfortran -c -o $@ $<

lapack_fortran/cwrappers/dgemm.o: lapack_fortran/cwrappers/kinds.o
lapack_fortran/cwrappers/sgemm.o: lapack_fortran/cwrappers/kinds.o

# C files
# lapack_c/cppwrappers/*.o: lapack_c/cppwrappers/*.cpp
# 	g++ -std=c++17 -c -o $@ $<

lapack_c/cppwrappers/gemm.o: lapack_c/cppwrappers/gemm.cpp
	g++ -std=c++17 -I./lapack_c/include -I./lapack_cpp/include -c -o $@ $<

# C++ files
lapack_cpp/cppwrappers/blas/*.o: lapack_cpp/cppwrappers/blas/*.cpp
	g++ -std=c++17 -c -o $@ $<

# Test files
example/test_gemm_fortran: example/test_gemm.f90 lapack_fortran/src/sgemm.o lapack_fortran/src/lsame.o  lapack_fortran/src/xerbla.o
	gfortran -o $@ $^

example/test_gemm_cpp: example/test_gemm.f90 lapack_fortran/cwrappers/sgemm.o lapack_c/cppwrappers/gemm.o
	gfortran -o $@ $^ -lstdc++

