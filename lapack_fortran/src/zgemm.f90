#ifndef USE_FORTRAN_BLAS
! Wrapper around zgemm defined in C
subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    use iso_c_binding, only: c_char
    use kinds, only: idx_kind
    implicit none
    ! Arguments
    double complex, intent(in) :: alpha,beta
    integer(idx_kind), intent(in) :: m,n,k,lda,ldb,ldc
    character, intent(in) :: transa,transb
    double complex, intent(in) :: a(lda,*),b(ldb,*)
    double complex, intent(inout) :: c(ldc,*)

    ! Interface to C
    interface
        subroutine lapack_c_zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(c, name="lapack_c_zgemm")
            use iso_c_binding, only: c_char, c_double_complex
            use kinds, only: idx_kind
            implicit none
            complex(c_double_complex), value, intent(in) :: alpha,beta
            integer(idx_kind), value, intent(in) :: m,n,k,lda,ldb,ldc
            character(c_char), value, intent(in) :: transa,transb
            complex(c_double_complex), intent(in) :: a(lda,*),b(ldb,*)
            complex(c_double_complex), intent(inout) :: c(ldc,*)
        end subroutine lapack_c_zgemm
    end interface

    character(c_char) :: transa_to_c, transb_to_c

    ! Convert characters to C characters
    transa_to_c = transa
    transb_to_c = transb

    ! Call C function
    call lapack_c_zgemm(transa_to_c, transb_to_c, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

end subroutine zgemm
#endif