! Wrapper around cgemm defined in C
subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    use iso_c_binding, only: c_char
    implicit none
    ! Arguments
    complex, intent(in) :: alpha,beta
    integer, intent(in) :: m,n,k,lda,ldb,ldc
    character, intent(in) :: transa,transb
    complex, intent(in) :: a(lda,*),b(ldb,*)
    complex, intent(inout) :: c(ldc,*)

    ! Interface to C
    interface
        subroutine lapack_c_cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(c, name="lapack_c_cgemm")
            use iso_c_binding, only: c_char, c_float_complex, c_int
            implicit none
            complex(c_float_complex), value, intent(in) :: alpha,beta
            integer(c_int), value, intent(in) :: m,n,k,lda,ldb,ldc
            character(c_char), value, intent(in) :: transa,transb
            complex(c_float_complex), intent(in) :: a(lda,*),b(ldb,*)
            complex(c_float_complex), intent(inout) :: c(ldc,*)
        end subroutine lapack_c_cgemm
    end interface

    character(c_char) :: transa_to_c, transb_to_c

    ! Convert characters to C characters
    transa_to_c = transa
    transb_to_c = transb

    ! Call C function
    call lapack_c_cgemm(transa_to_c, transb_to_c, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

end subroutine cgemm
