#include "lapack_c/gemm.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

//
// Define headers for fortran routines
//

#define FORTRAN_DGEMM FORTRAN_TO_C_NAME(dgemm, DGEMM)
void FORTRAN_DGEMM(const char* transa,
                   const char* transb,
                   const int* m,
                   const int* n,
                   const int* k,
                   const double* alpha,
                   const double* A,
                   const int* lda,
                   const double* B,
                   const int* ldb,
                   const double* beta,
                   double* C,
                   const int* ldc);

#define FORTRAN_SGEMM FORTRAN_TO_C_NAME(sgemm, SGEMM)
void FORTRAN_SGEMM(const char* transa,
                   const char* transb,
                   const int* m,
                   const int* n,
                   const int* k,
                   const float* alpha,
                   const float* A,
                   const int* lda,
                   const float* B,
                   const int* ldb,
                   const float* beta,
                   float* C,
                   const int* ldc);

#define FORTRAN_CGEMM FORTRAN_TO_C_NAME(cgemm, CGEMM)
void FORTRAN_CGEMM(const char* transa,
                   const char* transb,
                   const int* m,
                   const int* n,
                   const int* k,
                   const lapack_float_complex* alpha,
                   const lapack_float_complex* A,
                   const int* lda,
                   const lapack_float_complex* B,
                   const int* ldb,
                   const lapack_float_complex* beta,
                   lapack_float_complex* C,
                   const int* ldc);

#define FORTRAN_ZGEMM FORTRAN_TO_C_NAME(zgemm, ZGEMM)
void FORTRAN_ZGEMM(const char* transa,
                   const char* transb,
                   const int* m,
                   const int* n,
                   const int* k,
                   const lapack_double_complex* alpha,
                   const lapack_double_complex* A,
                   const int* lda,
                   const lapack_double_complex* B,
                   const int* ldb,
                   const lapack_double_complex* beta,
                   lapack_double_complex* C,
                   const int* ldc);

//
// Define actual wrappers
//
void lapack_c_dgemm(char transa,
               char transb,
               int m,
               int n,
               int k,
               double alpha,
               const double* A,
               int lda,
               const double* B,
               int ldb,
               double beta,
               double* C,
               int ldc)
{
    FORTRAN_DGEMM(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
                  C, &ldc);
}

void lapack_c_sgemm(char transa,
               char transb,
               int m,
               int n,
               int k,
               float alpha,
               const float* A,
               int lda,
               const float* B,
               int ldb,
               float beta,
               float* C,
               int ldc)
{
    FORTRAN_SGEMM(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
                  C, &ldc);
}

void lapack_c_cgemm(char transa,
                    char transb,
                    lapack_idx m,
                    lapack_idx n,
                    lapack_idx k,
                    lapack_float_complex alpha,
                    const lapack_float_complex* A,
                    lapack_idx lda,
                    const lapack_float_complex* B,
                    lapack_idx ldb,
                    lapack_float_complex beta,
                    lapack_float_complex* C,
                    lapack_idx ldc)
{
    FORTRAN_CGEMM(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
                  C, &ldc);
}

void lapack_c_zgemm(char transa,
                    char transb,
                    lapack_idx m,
                    lapack_idx n,
                    lapack_idx k,
                    lapack_double_complex alpha,
                    const lapack_double_complex* A,
                    lapack_idx lda,
                    const lapack_double_complex* B,
                    lapack_idx ldb,
                    lapack_double_complex beta,
                    lapack_double_complex* C,
                    lapack_idx ldc)
{
    FORTRAN_ZGEMM(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
                  C, &ldc);
}

#ifdef __cplusplus
}
#endif  // __cplusplus