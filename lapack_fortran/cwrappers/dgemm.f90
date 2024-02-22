! Wrapper around dgemm defined in C
subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    use iso_c_binding, only: c_char
    use kinds, only: idx_kind
    implicit none
    ! Arguments
    double precision, intent(in) :: alpha,beta
    integer(idx_kind), intent(in) :: m,n,k,lda,ldb,ldc
    character, intent(in) :: transa,transb
    double precision, intent(in) :: a(lda,*),b(ldb,*)
    double precision, intent(inout) :: c(ldc,*)

    ! Interface to C
    interface
        subroutine lapack_c_dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(c, name="lapack_c_dgemm")
            use iso_c_binding, only: c_char, c_double
            use kinds, only: idx_kind
            implicit none
            real(c_double), value, intent(in) :: alpha,beta
            integer(idx_kind), value, intent(in) :: m,n,k,lda,ldb,ldc
            character(c_char), value, intent(in) :: transa,transb
            real(c_double), intent(in) :: a(lda,*),b(ldb,*)
            real(c_double), intent(inout) :: c(ldc,*)
        end subroutine lapack_c_dgemm
    end interface

    character(c_char) :: transa_to_c, transb_to_c

    ! Convert characters to C characters
    transa_to_c = transa
    transb_to_c = transb

    ! Call C function
    call lapack_c_dgemm(transa_to_c, transb_to_c, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

end subroutine dgemm
