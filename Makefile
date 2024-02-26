CXX=g++
CXXFLAGS=-std=c++17 -Wall -Wextra -pedantic -g -DUSE_FORTRAN_BLAS

FC=gfortran
FFLAGS=-Wall -Wextra -pedantic -g -DUSE_FORTRAN_BLAS

all: ./example/test_gemm_fortran ./example/test_gemm_cpp ./example/test_gemv_cpp ./example/test_gemv_fortran ./example/test_trsv_cpp

HEADERS = $(wildcard lapack_c/include/*.h) \
		  $(wildcard lapack_c/include/*/*.h) \
		  $(wildcard lapack_c/include/*/*/*.h) \
		  $(wildcard lapack_cpp/include/*.hpp) \
		  $(wildcard lapack_cpp/include/*/*.hpp) \
		  $(wildcard lapack_cpp/include/*/*/*.hpp)

OBJFILES = lapack_fortran/src/sgemm.o \
		   lapack_fortran/src/dgemm.o \
		   lapack_fortran/src/cgemm.o \
		   lapack_fortran/src/zgemm.o \
		   lapack_c/src/gemm_c.o \
		   lapack_c/src/gemv_c.o \
		   lapack_c/src/rot_c.o \
		   lapack_c/src/trsv_c.o \
		   lapack_cpp/src/gemm_cpp.o \
		   lapack_cpp/src/gemv_cpp.o \
		   lapack_cpp/src/rot_cpp.o \
		   lapack_cpp/src/trsv_cpp.o \
		   lapack_fortran/src/zrot.o \
		   lapack_fortran/src/crot.o \

# Fortran files

# Note: The linker can't find zrot and crot unless they are recompiled.
lapack_fortran/src/crot.o: lapack_fortran/src/crot.f
	$(FC) $(FFLAGS) -cpp -c -o $@ $<

lapack_fortran/src/zrot.o: lapack_fortran/src/zrot.f
	$(FC) $(FFLAGS) -cpp -c -o $@ $<

lapack_fortran/src/%.o: lapack_fortran/src/%.f90
	$(FC) $(FFLAGS) -cpp -c -o $@ $<

lapack_fortran/src/sgemm.o: lapack_fortran/src/kinds.o

lapack_fortran/src/dgemm.o: lapack_fortran/src/kinds.o

lapack_fortran/src/cgemm.o: lapack_fortran/src/kinds.o

lapack_fortran/src/zgemm.o: lapack_fortran/src/kinds.o

# C files

lapack_c/src/%.o: lapack_c/src/%.cpp
	$(CXX) $(CXXFLAGS) -I./lapack_c/include -I./lapack_cpp/include -c -o $@ $<

lapack_c/src/%.o: lapack_c/src/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

# C++ files
lapack_cpp/src/%.o: lapack_cpp/src/%.cpp
	$(CXX) $(CXXFLAGS) -I./lapack_c/include -I./lapack_cpp/include -c -o $@ $<

# Test files

example/test_gemm_fortran: example/test_gemm.f90 $(OBJFILES) $(HEADERS)
	$(FC) $(FFLAGS) -o $@ $< $(OBJFILES) -lstdc++ -lblas

example/test_gemv_fortran: example/test_gemv.f90 $(OBJFILES) $(HEADERS)
	$(FC) $(FFLAGS) -o $@ $< $(OBJFILES) -lstdc++ -lblas

example/test_gemm_cpp: example/test_gemm.cpp $(OBJFILES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJFILES) -I./lapack_c/include -I./lapack_cpp/include -lstdc++ -lgfortran -lblas

example/test_gemv_cpp: example/test_gemv.cpp $(OBJFILES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJFILES) -I./lapack_c/include -I./lapack_cpp/include -lstdc++ -lgfortran -lblas

example/test_trsv_cpp: example/test_trsv.cpp $(OBJFILES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJFILES) -I./lapack_c/include -I./lapack_cpp/include -lstdc++ -lgfortran -lblas

clean:
	find . -type f -name '*.o' -delete