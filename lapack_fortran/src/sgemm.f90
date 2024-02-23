#ifndef USE_FORTRAN_BLAS
! Wrapper around sgemm defined in C
subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    use iso_c_binding, only: c_char
    use kinds, only: idx_kind
    implicit none
    ! Arguments
    real, intent(in) :: alpha,beta
    integer(idx_kind), intent(in) :: m,n,k,lda,ldb,ldc
    character, intent(in) :: transa,transb
    real, intent(in) :: a(lda,*),b(ldb,*)
    real, intent(inout) :: c(ldc,*)

    ! Interface to C
    interface
        subroutine lapack_c_sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(c, name="lapack_c_sgemm")
            use iso_c_binding, only: c_char, c_float
            use kinds, only: idx_kind
            implicit none
            real(c_float), value, intent(in) :: alpha,beta
            integer(idx_kind), value, intent(in) :: m,n,k,lda,ldb,ldc
            character(c_char), value, intent(in) :: transa,transb
            real(c_float), intent(in) :: a(lda,*),b(ldb,*)
            real(c_float), intent(inout) :: c(ldc,*)
        end subroutine lapack_c_sgemm
    end interface

    character(c_char) :: transa_to_c, transb_to_c

    ! Convert characters to C characters
    transa_to_c = transa
    transb_to_c = transb

    ! Call C function
    call lapack_c_sgemm(transa_to_c, transb_to_c, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

end subroutine sgemm
#endif