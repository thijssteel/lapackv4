! Define trsv on the fortran side using bind(c)
subroutine lapack_c_dtrsv(uplo, trans, diag, n, a, lda, x, incx) bind(c, name='lapack_c_dtrsv')
    use iso_c_binding
    implicit none
    character(kind=c_char), value, intent(in) :: uplo, trans, diag
    integer(kind=c_int), value, intent(in) :: n, lda, incx
    real(kind=c_double), dimension(lda, *), intent(in) :: a
    real(kind=c_double), dimension(*), intent(inout) :: x

    external dtrsv

    ! Forward call to fortran implementation
    call dtrsv( uplo, trans, diag, n, a, lda, x, incx )

end subroutine lapack_c_dtrsv

subroutine lapack_c_strsv(uplo, trans, diag, n, a, lda, x, incx) bind(c, name='lapack_c_strsv')
    use iso_c_binding
    implicit none
    character(kind=c_char), value, intent(in) :: uplo, trans, diag
    integer(kind=c_int), value, intent(in) :: n, lda, incx
    real(kind=c_float), dimension(lda, *), intent(in) :: a
    real(kind=c_float), dimension(*), intent(inout) :: x

    external strsv

    ! Forward call to fortran implementation
    call strsv( uplo, trans, diag, n, a, lda, x, incx )

end subroutine lapack_c_strsv

subroutine lapack_c_ctrsv(uplo, trans, diag, n, a, lda, x, incx) bind(c, name='lapack_c_ctrsv')
    use iso_c_binding
    implicit none
    character(kind=c_char), value, intent(in) :: uplo, trans, diag
    integer(kind=c_int), value, intent(in) :: n, lda, incx
    complex(kind=c_float_complex), dimension(lda, *), intent(in) :: a
    complex(kind=c_float_complex), dimension(*), intent(inout) :: x

    external ctrsv

    ! Forward call to fortran implementation
    call ctrsv( uplo, trans, diag, n, a, lda, x, incx )

end subroutine lapack_c_ctrsv

subroutine lapack_c_ztrsv(uplo, trans, diag, n, a, lda, x, incx) bind(c, name='lapack_c_ztrsv')
    use iso_c_binding
    implicit none
    character(kind=c_char), value, intent(in) :: uplo, trans, diag
    integer(kind=c_int), value, intent(in) :: n, lda, incx
    complex(kind=c_double_complex), dimension(lda, *), intent(in) :: a
    complex(kind=c_double_complex), dimension(*), intent(inout) :: x

    external ztrsv

    ! Forward call to fortran implementation
    call ztrsv( uplo, trans, diag, n, a, lda, x, incx )

end subroutine lapack_c_ztrsv

