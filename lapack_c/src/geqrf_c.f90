subroutine lapack_c_dgeqrf(m, n, A, lda, tau, work, lwork, info) bind(c, name='lapack_c_dgeqrf')
    use iso_c_binding
    implicit none
    integer(kind=c_int), value, intent(in) :: m, n, lda, lwork
    integer(kind=c_int), intent(out) :: info
    real(kind=c_double), dimension(*), intent(inout) :: a(lda, *), tau(*), work(*)

    external dgeqrf

    ! Forward call to fortran implementation
    call dgeqrf( m, n, A, lda, tau, work, lwork, info )

end subroutine lapack_c_dgeqrf

subroutine lapack_c_sgeqrf(m, n, A, lda, tau, work, lwork, info) bind(c, name='lapack_c_sgeqrf')
    use iso_c_binding
    implicit none
    integer(kind=c_int), value, intent(in) :: m, n, lda, lwork
    integer(kind=c_int), intent(out) :: info
    real(kind=c_float), dimension(*), intent(inout) :: a(lda, *), tau(*), work(*)

    external sgeqrf

    ! Forward call to fortran implementation
    call sgeqrf( m, n, A, lda, tau, work, lwork, info )

end subroutine lapack_c_sgeqrf

subroutine lapack_c_cgeqrf(m, n, A, lda, tau, work, lwork, info) bind(c, name='lapack_c_cgeqrf')
    use iso_c_binding
    implicit none
    integer(kind=c_int), value, intent(in) :: m, n, lda, lwork
    integer(kind=c_int), intent(out) :: info
    real(kind=c_float_complex), dimension(*), intent(inout) :: a(lda, *), tau(*), work(*)

    external cgeqrf

    ! Forward call to fortran implementation
    call cgeqrf( m, n, A, lda, tau, work, lwork, info )

end subroutine lapack_c_cgeqrf

subroutine lapack_c_zgeqrf(m, n, A, lda, tau, work, lwork, info) bind(c, name='lapack_c_zgeqrf')
    use iso_c_binding
    implicit none
    integer(kind=c_int), value, intent(in) :: m, n, lda, lwork
    integer(kind=c_int), intent(out) :: info
    real(kind=c_double_complex), dimension(*), intent(inout) :: a(lda, *), tau(*), work(*)

    external zgeqrf

    ! Forward call to fortran implementation
    call zgeqrf( m, n, A, lda, tau, work, lwork, info )

end subroutine lapack_c_zgeqrf